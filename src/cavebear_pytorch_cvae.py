## ========================================
## new version modified based on original version
## modified px_r to be gene and cell specific, modified #layer and #hidden
"""
Conditional VAE for scRNA-seq with ZINB loss (scVI-like)
Includes:
 - Species + Batch conditioning
 - ZINB loss with library size scaling
 - Sigmoid KL warmup schedule
 - Early stopping
 - LR scheduler (ReduceLROnPlateau)
 - Exports latent embeddings
 - Optional adversarial discriminator to align species in latent space
"""
import scanpy as sc
import numpy as np
import anndata
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, random_split, TensorDataset
from scipy import sparse
import pandas as pd
from pathlib import Path
import sys
import os, argparse, json
from sklearn.neighbors import NearestNeighbors
import random
import math
import matplotlib.pyplot as plt
import json

# ---------- Utilities ----------
def ensure_sparse_csr(X):
    """Ensure X is a scipy.sparse CSR matrix"""
    if sparse.isspmatrix_csr(X):
        return X
    elif sparse.issparse(X):
        return X.tocsr()
    else:
        return sparse.csr_matrix(X)

def build_dataset(adata, species_col, batch_col, dis):
    """
    Shared helper to build cond_matrix, size_factors and AnnDataset.
    Returns: dataset, species_cats, species_idx, cond_dim
    """
    libsize = np.asarray(adata.X.sum(axis=1)).squeeze().astype(np.float32)
    libsize[libsize == 0] = 1.0

    species = adata.obs[species_col].astype(str).values
    batch   = adata.obs[batch_col].astype(str).values

    species_cats, species_idx = np.unique(species, return_inverse=True)
    batch_cats,   batch_idx   = np.unique(batch,   return_inverse=True)

    cond_matrix = np.concatenate(
        [np.eye(len(species_cats))[species_idx],
         np.eye(len(batch_cats))[batch_idx]], axis=1
    ).astype(np.float32)

    dataset = AnnDataset(adata, cond_matrix, libsize,
                         species_idx=species_idx if dis else None)

    return dataset, species_cats, species_idx, cond_matrix.shape[1]


# ---------- Dataset ----------
class AnnDataset(Dataset):
    def __init__(self, adata, cond_matrix, size_factors, species_idx=None):
        # Keep sparse representation
        self.X = adata.X
        
        # Other dense arrays
        self.cond = torch.from_numpy(cond_matrix.astype(np.float32))
        self.size_factors = torch.from_numpy(size_factors.astype(np.float32)).unsqueeze(1)
        # species_idx: integer label per cell for discriminator
        self.species_idx = torch.from_numpy(species_idx.astype(np.int64)) if species_idx is not None else None


    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        # Convert one sparse row to dense float32 tensor
        x = self.X[idx].toarray().squeeze()
        x = torch.from_numpy(x)
        out = {"x": x, "cond": self.cond[idx], "sf": self.size_factors[idx]}
        if self.species_idx is not None:
            out["species_idx"] = self.species_idx[idx]
        return out

# ---------- Model ----------
class Encoder(nn.Module):
    def __init__(self, n_input, cond_dim, hidden_size, latent_dim, n_layers=3, dropout=0.1, negative_slope=0.2):
        super().__init__()
        layers = []
        input_dim = n_input + cond_dim
        for i in range(n_layers):
            layers.append(nn.Linear(input_dim if i == 0 else hidden_size, hidden_size))
            layers.append(nn.LayerNorm(hidden_size, elementwise_affine=True))
            layers.append(nn.LeakyReLU(negative_slope))
            layers.append(nn.Dropout(dropout))
        self.net = nn.Sequential(*layers)
        self.fc_mu = nn.Linear(hidden_size, latent_dim)
        self.fc_logvar = nn.Linear(hidden_size, latent_dim)

    def forward(self, x, c):
        x = x.float() # added to avoid RunTimeError: mat1 and mat2 must have the same dtype, but got Double and Float
        c = c.float() # added to avoid RunTimeError: mat1 and mat2 must have the same dtype, but got Double and Float
        h = torch.cat([x, c], dim=-1)
        h = self.net(h)
        return self.fc_mu(h), self.fc_logvar(h)

class Decoder(nn.Module):
    def __init__(self, latent_dim, cond_dim, hidden_size, n_genes, n_layers=3, dropout=0.1, negative_slope=0.2):
        super().__init__()
        layers = []
        input_dim = latent_dim + cond_dim
        for i in range(n_layers):
            layers.append(nn.Linear(input_dim if i == 0 else hidden_size, hidden_size))
            layers.append(nn.LayerNorm(hidden_size, elementwise_affine=True))
            layers.append(nn.LeakyReLU(negative_slope))
            layers.append(nn.Dropout(dropout))
        self.net = nn.Sequential(*layers)

        self.px_scale   = nn.Linear(hidden_size, n_genes)    # logits before softmax
        self.px_dropout = nn.Linear(hidden_size, n_genes)  # logits for dropout (pi)
        self.px_r_lin   = nn.Linear(hidden_size, n_genes)    # raw dispersion -> make positive in forward

        self.softmax  = nn.Softmax(dim=-1)
        self.softplus = nn.Softplus()

    def forward(self, z, c):
        h = torch.cat([z, c], dim=-1)
        h = self.net(h)
        px_decoder = self.softmax(self.px_scale(h))
        px_dropout = self.px_dropout(h)
        theta      = self.softplus(self.px_r_lin(h)) + 1e-8
        return px_decoder, theta, px_dropout

class Discriminator(nn.Module):
    def __init__(self, latent_dim, n_species, dropout=0.1):
        super().__init__()
        hidden = max(int((latent_dim * n_species) ** 0.5), 16)
        self.net = nn.Sequential(
            nn.Linear(latent_dim, hidden),
            nn.LayerNorm(hidden),
            nn.LeakyReLU(0.2),
            nn.Dropout(dropout),
            nn.Linear(hidden, n_species),
        )

    def forward(self, z):
        return self.net(z)

class CVAE(nn.Module):
    def __init__(self, n_genes, cond_dim, hidden_size, latent_dim,
                 n_layers=3, n_species=0, discriminator_weight=1.0):
        super().__init__()
        self.encoder = Encoder(n_genes, cond_dim, hidden_size, latent_dim, n_layers)
        self.decoder = Decoder(latent_dim, cond_dim, hidden_size, n_genes, n_layers)
        self.discriminator_weight = discriminator_weight
        self.discriminator = Discriminator(latent_dim, n_species) if n_species > 1 else None

    def reparam(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x_counts, cond, size_factors):
        x_enc = torch.log1p(x_counts)
        q_mu, q_logvar = self.encoder(x_enc, cond)
        z = self.reparam(q_mu, q_logvar)
        px_decoder, theta, px_dropout = self.decoder(z, cond)
        mu = px_decoder * size_factors
        pi = torch.sigmoid(px_dropout)
        return {"mu": mu, "theta": theta, "pi": pi,
                "q_mu": q_mu, "q_logvar": q_logvar, "z": z}


# ---------- ZINB & KL ----------
def log_zinb_positive(x, mu, theta, eps=1e-8):
    if theta.ndim == 1:
        theta = theta.unsqueeze(0)
    t1 = torch.lgamma(theta + x + 1e-8) - torch.lgamma(theta + 1e-8) - torch.lgamma(x + 1.0)
    t2 = theta * (torch.log(theta + eps) - torch.log(mu + theta + eps))
    t3 = x * (torch.log(mu + eps) - torch.log(mu + theta + eps))
    return t1 + t2 + t3


def zinb_loss(x, mu, theta, pi, eps=1e-8):
    log_nb = log_zinb_positive(x, mu, theta, eps)
    nb_case_zero    = torch.exp(log_nb)
    log_prob_zero   = torch.log(pi + (1.0 - pi) * nb_case_zero + eps)
    log_prob_nonzero = torch.log(1.0 - pi + eps) + log_nb
    mask_zero = (x <= 0.0).float()
    log_prob  = mask_zero * log_prob_zero + (1.0 - mask_zero) * log_prob_nonzero
    neg_log_likelihood = -torch.sum(log_prob, dim=-1)
    return neg_log_likelihood.mean(), neg_log_likelihood


def kl_divergence_gauss(q_mu, q_logvar):
    return -0.5 * torch.sum(1 + q_logvar - q_mu.pow(2) - q_logvar.exp(), dim=1)


# ---------- KL warmup (linear schedule) ----------
def kl_weight_schedule(epoch, steepness=0.0025):
    return float(steepness * epoch)


# ---------- Training helpers ----------
def compute_loss(model, batch, kl_weight, device, species_idx=None):
    x    = batch["x"].to(device).float()
    cond = batch["cond"].to(device).float()
    sf   = batch["sf"].to(device).float()
    out  = model(x, cond, sf)

    mu, theta, pi   = out["mu"], out["theta"], out["pi"]
    q_mu, q_logvar  = out["q_mu"], out["q_logvar"]

    recon_loss, _  = zinb_loss(x, mu, theta, pi)
    kl             = kl_divergence_gauss(q_mu, q_logvar).mean()
    total_loss     = recon_loss + kl_weight * kl

    disc_loss = torch.tensor(0.0, device=device)
    if model.discriminator is not None and species_idx is not None:
        logits    = model.discriminator(out["q_mu"])
        disc_loss = F.cross_entropy(logits, species_idx.to(device))

    return total_loss, recon_loss, kl, disc_loss

def evaluate(model, dataloader, device, kl_weight):
    """Evaluate on a dataloader; returns (total, recon, kl, disc) losses."""
    model.eval()
    total_loss = total_recon = total_kl = total_disc = 0.0
    with torch.no_grad():
        for batch in dataloader:
            species_idx = batch.get("species_idx")
            loss, recon, kl, disc_loss = compute_loss(
                model, batch, kl_weight, device, species_idx)
            n = batch["x"].size(0)
            total_loss  += loss.item()      * n
            total_recon += recon.item()     * n
            total_kl    += kl.item()        * n
            total_disc  += disc_loss.item() * n
    n = len(dataloader.dataset)
    return total_loss / n, total_recon / n, total_kl / n, total_disc / n


def collate_dense(batch):
    out = {
        "x":    torch.stack([b["x"]    for b in batch]),
        "cond": torch.stack([b["cond"] for b in batch]),
        "sf":   torch.stack([b["sf"]   for b in batch]),
    }
    if "species_idx" in batch[0]:
        out["species_idx"] = torch.stack([b["species_idx"] for b in batch])
    return out


def setTrainingSets(dataset, batch_size, seed):
    random.seed(seed)
    val_split = 0.1
    val_size   = int(val_split * len(dataset))
    train_size = len(dataset) - val_size
    print(f"Validation Set Size: {val_size}, Training Set Size: {train_size}")

    train_set, val_set = random_split(dataset, [train_size, val_size])

    train_loader = DataLoader(train_set, batch_size=batch_size, shuffle=True,  collate_fn=collate_dense)
    val_loader   = DataLoader(val_set,   batch_size=batch_size, shuffle=False, collate_fn=collate_dense)
    return train_loader, val_loader, val_set, train_size

# ---------- fit_CVAE ----------
def fit_CVAE(model, lr, patience, train_loader, val_loader, train_size,
             outdir, model_file, device, epochs=500,
             dis=False, discriminator_weight=1.0):

    # VAE params only (encoder + decoder); discriminator has its own optimizer
    vae_params = list(model.encoder.parameters()) + list(model.decoder.parameters())
    optimizer  = torch.optim.Adam(vae_params, lr=lr)
    scheduler  = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=patience, min_lr=1e-5)

    optimizer_disc = None
    if dis and model.discriminator is not None:
        optimizer_disc = torch.optim.Adam(model.discriminator.parameters(), lr=0.0001)

    best_val_loss   = float("inf")
    patience_counter = 0

    for epoch in range(1, epochs + 1):
        kl_weight = kl_weight_schedule(epoch)
        model.train()
        total_loss = total_recon = total_kl = 0.0

        for batch in train_loader:
            species_idx = batch.get("species_idx")

            # --- Step 1: update discriminator (maximise classification accuracy) ---
            if optimizer_disc is not None and species_idx is not None:
                optimizer_disc.zero_grad()
                x    = batch["x"].to(device).float()
                cond = batch["cond"].to(device).float()
                sf   = batch["sf"].to(device).float()
                with torch.no_grad():
                    out = model(x, cond, sf)
                logits    = model.discriminator(out["q_mu"].detach())
                disc_loss = F.cross_entropy(logits, species_idx.to(device))
                disc_loss.backward()
                optimizer_disc.step()

            # --- Step 2: update VAE (minimise recon+KL, maximise disc confusion) ---
            optimizer.zero_grad()
            loss, recon, kl, disc_loss = compute_loss(
                model, batch, kl_weight, device, species_idx)
            if dis and model.discriminator is not None:
                loss = loss - discriminator_weight * disc_loss
            loss.backward()
            optimizer.step()

            n = batch["x"].size(0)
            total_loss  += loss.item()  * n
            total_recon += recon.item() * n
            total_kl    += kl.item()    * n

        train_loss = total_loss / train_size

        val_loss, val_recon, val_kl, val_disc = evaluate(model, val_loader, device, kl_weight)
        val_score = val_recon - discriminator_weight * val_disc if dis else val_recon

        print(
            f"Epoch {epoch:03d} | KLw {kl_weight:.3f} | "
            f"Train {train_loss:.4f} | Val {val_loss:.4f} | "
            f"Recon {val_recon:.4f} | KL {val_kl:.4f} | "
            f"Disc {val_disc:.4f} | LR {optimizer.param_groups[0]['lr']:.2e}"
        )

        scheduler.step(val_score)

        if val_score < best_val_loss:
            best_val_loss    = val_score
            patience_counter = 0
            torch.save(model.state_dict(), str(outdir) + '/' + str(model_file))
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch} "
                      f"(no improvement for {patience} epochs)")
                break

    return model, best_val_loss

# ---------- LISI ----------
def compute_lisi(X, label, perplexity=30):
    """
    Compute the Local Inverse Simpson Index (LISI).
    Adapted from https://github.com/slowkow/harmonypy/blob/master/harmonypy/lisi.py
    """
    knn = NearestNeighbors(n_neighbors=perplexity * 3, algorithm='kd_tree').fit(X)
    distances, indices = knn.kneighbors(X)
    indices   = indices[:, 1:]
    distances = distances[:, 1:]
    labels       = pd.Categorical(label)
    n_categories = len(labels.categories)
    print(" Computing simpson val:")
    simpson = compute_simpson(distances.T, indices.T, labels, n_categories, perplexity)
    return np.mean(1 / simpson)


def compute_simpson(distances, indices, labels, n_categories, perplexity, tol=1e-5):
    """Adapted from harmonypy."""
    n       = distances.shape[1]
    P       = np.zeros(distances.shape[0])
    simpson = np.zeros(n)
    logU    = np.log(perplexity)

    for i in range(n):
        beta    = 1
        betamin = -np.inf
        betamax =  np.inf

        P     = np.exp(-distances[:, i] * beta)
        P_sum = np.sum(P)
        if P_sum == 0:
            H = 0
            P = np.zeros(distances.shape[0])
        else:
            H = np.log(P_sum) + beta * np.sum(distances[:, i] * P) / P_sum
            P = P / P_sum
        Hdiff = H - logU

        for t in range(50):
            if abs(Hdiff) < tol:
                break
            if Hdiff > 0:
                betamin = beta
                beta = beta * 2 if not np.isfinite(betamax) else (beta + betamax) / 2
            else:
                betamax = beta
                beta = beta / 2 if not np.isfinite(betamin) else (beta + betamin) / 2

            P     = np.exp(-distances[:, i] * beta)
            P_sum = np.sum(P)
            if P_sum == 0:
                H = 0
                P = np.zeros(distances.shape[0])
            else:
                H = np.log(P_sum) + beta * np.sum(distances[:, i] * P) / P_sum
                P = P / P_sum
            Hdiff = H - logU

        if H == 0:
            simpson[i] = -1
        for label_category in labels.categories:
            ix    = indices[:, i]
            q     = labels[ix] == label_category
            if np.any(q):
                P_sum = np.sum(P[q])
                simpson[i] += P_sum * P_sum

    return simpson

# ---------- Time Predictor ----------
class TimeMLP(nn.Module):
    def __init__(self, input_dim, embed_dim, nlayer, dropout=0.1):
        super().__init__()
        assert nlayer >= 1
        layers = []
        for i in range(nlayer):
            layers.append(nn.Linear(input_dim if i == 0 else embed_dim, embed_dim))
            layers.append(nn.LayerNorm(embed_dim, elementwise_affine=True))
            layers.append(nn.LeakyReLU(negative_slope=0.01, inplace=True))
            layers.append(nn.Dropout(dropout))
        self.backbone = nn.Sequential(*layers)
        self.out = nn.Linear(embed_dim, 1)

    def forward(self, x):
        return F.relu(self.out(self.backbone(x)))


class TimePredictor:
    def __init__(self, input_dim, embed_dim, nlayer, output_model,
                 learning_rate_x, device=None):
        self.learning_rate_x = learning_rate_x
        self.nlayer      = nlayer
        self.input_dim   = input_dim
        self.embed_dim   = embed_dim
        self.output_model = output_model

        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(device)

        self.model = TimeMLP(input_dim=self.input_dim,
                             embed_dim=self.embed_dim,
                             nlayer=self.nlayer).to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(),
                                          lr=self.learning_rate_x, eps=1e-5)

    def _loss(self, y_true, y_pred):
        return F.mse_loss(y_pred.view(-1), y_true.view(-1))

    @torch.no_grad()
    def get_losses_rna(self, data_x, label_x):
        self.model.eval()
        x    = torch.as_tensor(data_x,  dtype=torch.float32, device=self.device)
        y    = torch.as_tensor(label_x, dtype=torch.float32, device=self.device)
        pred = self.model(x).squeeze(1)
        return float(self._loss(y, pred).detach().cpu().item())

    @torch.no_grad()
    def predict_time(self, data_x):
        self.model.eval()
        x = torch.as_tensor(data_x, dtype=torch.float32, device=self.device)
        return self.model(x).detach().cpu().numpy()

    def save(self, folder):
        os.makedirs(folder, exist_ok=True)
        torch.save({
            "state_dict": self.model.state_dict(),
            "optimizer":  self.optimizer.state_dict(),
            "config": {
                "input_dim":      self.input_dim,
                "embed_dim":      self.embed_dim,
                "nlayer":         self.nlayer,
                "learning_rate_x": self.learning_rate_x,
            },
        }, os.path.join(folder, "mymodel.pt"))

    def restore(self, restore_folder, load_optimizer=True):
        ckpt = torch.load(os.path.join(restore_folder, "mymodel.pt"),
                          map_location=self.device)
        self.model.load_state_dict(ckpt["state_dict"])
        if load_optimizer and "optimizer" in ckpt:
            self.optimizer.load_state_dict(ckpt["optimizer"])

    def train(self, data_x, data_x_val, label_x, label_x_val,
              output_model, batch_size, my_epochs, patience):
        if output_model:
            self.output_model = output_model

        iter_list       = []
        train_loss_list = []
        val_loss_list   = []

        ckpt_path = os.path.join(self.output_model, "mymodel.pt")
        if os.path.exists(ckpt_path):
            self.restore(self.output_model)
        else:
            loss_val_best = math.inf
            
            train_ds = TensorDataset(
                torch.as_tensor(data_x,     dtype=torch.float32),
                torch.as_tensor(label_x,    dtype=torch.float32))
            val_ds = TensorDataset(
                torch.as_tensor(data_x_val,  dtype=torch.float32),
                torch.as_tensor(label_x_val, dtype=torch.float32))

            train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True,  drop_last=False)
            val_loader   = DataLoader(val_ds,   batch_size=batch_size, shuffle=False, drop_last=False)

            last_improvement = 0

            for epoch in range(1, my_epochs):
                iter_list.append(epoch)
                self.model.train()

                for xb, yb in train_loader:
                    xb   = xb.to(self.device, non_blocking=True)
                    yb   = yb.to(self.device, non_blocking=True)
                    pred = self.model(xb).squeeze(1)
                    loss = self._loss(yb, pred)
                    self.optimizer.zero_grad(set_to_none=True)
                    loss.backward()
                    self.optimizer.step()

                self.model.eval()
                with torch.no_grad():
                    tlosses = [self._loss(
                        yb.to(self.device), self.model(xb.to(self.device)).squeeze(1)
                    ).item() for xb, yb in train_loader]
                    train_loss = float(np.nanmean(tlosses))
                    train_loss_list.append(train_loss)

                    vlosses = [self._loss(
                        yb.to(self.device), self.model(xb.to(self.device)).squeeze(1)
                    ).item() for xb, yb in val_loader]
                    val_loss = float(np.nanmean(vlosses))
                    val_loss_list.append(val_loss)

                print(f"epoch {epoch} | train_loss={train_loss:.6f} | val_loss={val_loss:.6f}")

                if np.isnan(val_loss):
                    break

                if val_loss <= loss_val_best:
                    self.save(self.output_model)
                    loss_val_best    = val_loss
                    last_improvement = 0
                else:
                    last_improvement += 1

                if last_improvement > patience:
                    if os.path.exists(ckpt_path):
                        self.restore(self.output_model, load_optimizer=False)
                    break

        return iter_list, train_loss_list, val_loss_list


def split_time_trainval(adata, train_species, seed, nsubsample=10000):
    """Split train/val/test cells for temporal prediction."""
    adata_train = adata[(adata.obs.species == train_species) &
                        (adata.obs.time_label == adata.obs.time_label)]
    X_train = adata_train.obsm["X_cvae"]

    random.seed(seed)
    adata_test_index = random.sample(
        range(adata_train.shape[0]),
        min(int(adata_train.shape[0] * 0.1), nsubsample))
    data_test       = X_train[adata_test_index, :]
    data_test_label = adata_train[adata_test_index].obs.time_label

    list_subset = list(set(range(adata_train.shape[0])) - set(adata_test_index))
    adata_train = adata_train[list_subset]
    X_train     = X_train[list_subset]

    adata_val_index = random.sample(
        range(adata_train.shape[0]),
        min(int(adata_train.shape[0] * 0.1), nsubsample))
    data_val       = X_train[adata_val_index, :]
    data_val_label = adata_train[adata_val_index].obs.time_label

    list_subset     = list(set(range(adata_train.shape[0])) - set(adata_val_index))
    data_train      = X_train[list_subset]
    data_train_label = adata_train[list_subset].obs.time_label

    return data_train, data_train_label, data_val, data_val_label, data_test, data_test_label

# ---------- Time prediction helper (for reuse without retraining) ----------
def predict_time_from_saved(time_model_outdir, adata_target, train_species, device=None):
    """
    Load the best saved temporal model and predict pseudotime for new cells.

    Parameters
    ----------
    time_model_outdir : str or Path
        Path to the time_prediction_models/{model_name}/ folder.
    adata_target : AnnData
        Cells to predict time for. Must have obsm['X_cvae'] already set.
    train_species : str
        Training species used when the model was saved (used to find best_model record).
    device : str or None
        'cuda', 'cpu', or None (auto-detect).

    Returns
    -------
    np.ndarray of shape (n_cells, 1) with predicted pseudotime values.

    Example
    -------
    preds = predict_time_from_saved(
        time_model_outdir='/path/to/time_prediction_models/cvae_pytorch_best_model_...',
        adata_target=adata_human,
        train_species='mouse')
    adata_human.obs['pred_time'] = preds
    """
    time_model_outdir = Path(time_model_outdir)
    best_model_record = time_model_outdir / f'best_model_{train_species}.json'

    if not best_model_record.exists():
        raise FileNotFoundError(
            f"Best model record not found: {best_model_record}\n"
            f"Run with predict='time' first to train and save the temporal model.")

    with open(best_model_record) as f:
        cfg = json.load(f)

    print(f"Loading best temporal model from {time_model_outdir / 'best_model'}")
    print(f"  Config: lr={cfg['lr']}  nlayer={cfg['nlayer']}  "
          f"embed_dim={cfg['embed_dim']}  input_dim={cfg['input_dim']}")

    predictor = TimePredictor(
        input_dim=cfg['input_dim'],
        embed_dim=cfg['embed_dim'],
        nlayer=cfg['nlayer'],
        output_model=str(time_model_outdir / 'best_model'),
        learning_rate_x=cfg['lr'],
        device=device)

    predictor.restore(str(time_model_outdir / 'best_model'), load_optimizer=False)
    return predictor.predict_time(adata_target.obsm["X_cvae"])

# ---------- parse hyperparameter string for LISI_log file ----------
def parse_hyperparameters(hyperparam_str):
    parts = hyperparam_str.split("_")

    parsed = {
        "lr": float(parts[0]),
        "n_layers": int(parts[1]),
        "latent_dim": int(parts[2]),
        "dis": 0.0,  # default if not present
        "seed": 101,  # default if not present
    }

    remaining = parts[3:]
    for part in remaining:
        if part.startswith("dis"):
            weight_str = part[3:]
            parsed["dis"] = float(weight_str) if weight_str else 0.0
        elif part.startswith("seed"):
            seed_str = part[4:]
            parsed["seed"] = int(seed_str) if seed_str else 101
        else:
            parsed["target_time"] = part

    return parsed


# ---------- Main ----------
def main(args):
    input_h5ad     = args.input_h5ad
    target_species = args.target_species
    train_species  = args.train_species
    target_time    = args.target_time
    group          = args.group
    batch_col      = args.batch_colname
    learning_rate  = args.learning_rate
    embed_dim      = args.embed_dim
    nlayer         = args.nlayer
    batch_size     = args.batch_size
    patience       = args.patience
    my_epochs      = args.my_epochs
    hidden_size    = args.hidden_size
    predict        = args.predict
    dis            = args.dis == 'dis'
    discriminator_weight = args.discriminator_weight
    time_label     = args.time_label
    seed           = args.seed

    ## hyperparameters
    lr        = float(learning_rate)
    n_layers  = int(nlayer)
    latent_dim = int(embed_dim)
    epochs    = my_epochs

    torch.set_default_dtype(torch.float32)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    ## set model parameters
    input_name = str(input_h5ad.split("/")[-1].split(".h5ad")[0])
    species_col = "species"

    # ── model name ────────────────────────────────────────────────────────
    hyperparameters = f"{str(lr)}_{str(n_layers)}_{str(latent_dim)}" # main hyperparameter settings
    if dis: # add discriminator to hyperparameters if it is used
        if discriminator_weight != 0.0:
            hyperparameters += f"_dis{discriminator_weight}"   # e.g. _dis2.0
    if target_time != '': # add target_time to hyperparameters if it is specified
        hyperparamters += f"_{target_time}"
    if seed != 101:
        hyperparameters += f"_seed{seed}"

    model_name = f"cvae_pytorch_disc_best_model_{str(input_name)}_{hyperparameters}"
    print(f"The model name is: {model_name}")

    model_file = str(model_name) + ".pth"
    results_dir = Path.cwd().parent / "results"
    results_dir.mkdir(parents=True, exist_ok=True) # make the results folder if it doesn't already exist

    # Check if the model already exists and if it does not, create the folder, else set outdir to the folder it exists in
    model_path = Path(f"{results_dir}/{input_name}/{hyperparameters}/")
    matches = list(model_path.glob(model_file))
    if len(matches) == 0:
        outdir = model_path
        os.makedirs(outdir, exist_ok=True)
        print("The model could not be found, training will be run.")
        print(f"Output Folder: {str(outdir)}")
    else:
        outdir = matches[0].parent
        print("The model was found in folder: " + str(outdir))

    ## check whether GPU is used
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Load the data and filter/ add obs labels as determined by arguments
    print("Loading AnnData:", input_h5ad)
    adata = anndata.read_h5ad(input_h5ad)
    print("AnnData was loaded!")
    if target_species != 'human':
        adata = adata[adata.obs['species'] != 'human']
    else:
        if 'origin' in adata.obs.columns:
            adata = adata[(adata.obs['species'] != 'zebrafish') &
                          (adata.obs['origin']  != 'Disteche')]
        else:
            adata = adata[(adata.obs['species'] != 'zebrafish')]

    adata.obs = adata.obs.copy()
    adata.obs['group'] = adata.obs[group] if (group != '' and group in adata.obs.columns) else ''
    adata.obs['time_label'] = adata.obs[time_label] if (time_label != '' and time_label in adata.obs.columns) else ''

    if target_time != '':
        adata = adata[(adata.obs.species != target_species) |
                      (adata.obs.time_label == float(target_time))]
    print(f'adata shape: {adata.shape}')

    embeddings_path = str(outdir) + "/" + model_name + "_cell_embeddings.npy" # set the cell embeddings file name

    if predict == 'train':
        if not Path(embeddings_path).exists():
            print(f"Training and/or evaluation of {model_name} is being performed!")

            dataset, species_cats, _, cond_dim = build_dataset(
                adata, species_col, batch_col, dis)

            n_cells, n_genes = adata.X.shape
            n_species = len(species_cats) if dis else 0
            print(f"Cells: {n_cells}, Genes: {n_genes}, Cond dim: {cond_dim}")
            print(f"Species: {species_cats} -> n_species for discriminator: {n_species}")

            train_loader, val_loader, data_val, train_size = setTrainingSets(dataset, batch_size, seed)

            if Path(str(outdir) + '/' + model_file).exists():
                print("The model already exists. Skipping training!")
                model = CVAE(n_genes, cond_dim, hidden_size, latent_dim, n_layers,
                            n_species=n_species,
                            discriminator_weight=discriminator_weight).to(device).float()
            else:
                print("The model does not exist. Training is being performed!")
                model = CVAE(n_genes, cond_dim, hidden_size, latent_dim, n_layers,
                            n_species=n_species,
                            discriminator_weight=discriminator_weight).to(device).float()
                model, best_val_loss = fit_CVAE(
                    model, lr, patience, train_loader, val_loader, train_size,
                    outdir, model_file, device, epochs,
                    dis=dis, discriminator_weight=discriminator_weight)
                print("Best validation loss:", best_val_loss)

            model.load_state_dict(torch.load(str(outdir) + '/' + model_file,
                                            map_location=device))
            model.eval()
            full_loader = DataLoader(dataset, batch_size=batch_size, shuffle=False,
                                    collate_fn=collate_dense)
            all_mu = []
            with torch.no_grad():
                for batch in full_loader:
                    x    = batch["x"].to(device).float()
                    cond = batch["cond"].to(device).float()
                    sf   = batch["sf"].to(device).float()
                    out  = model(x, cond, sf)
                    all_mu.append(out["q_mu"].cpu().numpy())

            latent_means = np.vstack(all_mu)
            np.save(embeddings_path, latent_means)
            print("Saved cell_embeddings.npy")

        else:
            print('Training was already performed. Loading model and computing UMAP/LISI evaluation!')
            latent_means = np.load(embeddings_path)

            dataset, species_cats, _, cond_dim = build_dataset(
                adata, species_col, batch_col, dis)
            _, _, data_val, _ = setTrainingSets(dataset, batch_size, seed)

        adata.obsm["X_cvae"] = latent_means

        # ── UMAP ──────────────────────────────────────────────────────────────
        umap_path = f'{outdir}/{str(model_name)}_umap.npy'
        if Path(umap_path).exists():
            print(f'Skipping UMAP Calculation. {umap_path} already exists.')
        else:
            print("Computing UMAP")
            # Make the UMAPs using scanpy
            sc.settings.figdir = f'{outdir}/umaps/'
            sc.settings.file_format_figs = "png"    # change extension globally
            # First create a UMAP for before CVAE training
            sc.pp.neighbors(adata, n_neighbors = 10)
            sc.tl.umap(adata)
            sc.pl.umap(adata, color=['species'], save=f'_{str(model_name)}_beforeTraining_species.png')
            # Second, create a UMAP for after CVAE training
            sc.pp.neighbors(adata, n_neighbors = 10, use_rep= 'X_cvae')
            sc.tl.umap(adata)
            np.save(f'{outdir}/{str(model_name)}_umap.npy', adata.obsm['X_umap'])
            sc.pl.umap(adata, color=['species'], save=f'_{str(model_name)}_species.png')
            sc.pl.umap(adata, color=['batch'], save=f'_{str(model_name)}_batch.png')
            sc.pl.umap(adata, color=['group'], save=f'_{str(model_name)}_group.png')
            # plot species individually
            species_list = adata.obs['species'].unique()
            n_species = len(species_list)
            colors = sc.pl.palettes.default_20[:n_species]
            colors[0], colors[1] = colors[1], colors[0]
            palette = dict(zip(adata.obs['species'].cat.categories, colors))
            fig, axes = plt.subplots(1, n_species, figsize=(6 * n_species, 5))
            if n_species == 1:
                axes = [axes]
            for ax, sp in zip(axes, species_list):
                sc.pl.umap(
                    adata,
                    color='species',
                    groups=[sp],       # only highlight this species
                    palette=palette,
                    ax=ax,
                    size=1,
                    show=False,
                    title=sp,
                    na_color='white', 
                    legend_loc=None
                )
                for collection in ax.collections:
                    collection.set_rasterized(True)
            plt.tight_layout()
            fig.savefig(f"{outdir}/umaps/umap_{str(model_name)}_species_sep.png", dpi=300, bbox_inches='tight', transparent=True)
            fig.savefig(f"{outdir}/umaps/umap_{str(model_name)}_species_sep.pdf", dpi=300, bbox_inches='tight', transparent=True)

        # ── LISI ──────────────────────────────────────────────────────────────
        print("Computing LISI")
        idx        = data_val.indices
        val_subset = adata[idx].copy()
        X_val      = np.asarray(val_subset.obsm["X_cvae"])
        # Standardize each dimension of embedding
        X_val      = (X_val - X_val.mean(axis=0)) / (X_val.std(axis=0))
        val_lisi   = compute_lisi(X_val, val_subset.obs.species.tolist(), perplexity=30)
        
        # Add LISI score for hyperparameter settings to the log
        parsed_params = parse_hyperparameters(hyperparameters)
        new_row = pd.DataFrame([{**parsed_params, "lisi_score": val_lisi}])
        lisi_log = Path(outdir).parent / f"LISI_log.txt"
        if os.path.exists(lisi_log):
            existing_df = pd.read_csv(lisi_log)
            hyperparam_cols = [col for col in existing_df.columns if col != "lisi_score"]
            
            duplicate_mask = (existing_df[hyperparam_cols] == new_row[hyperparam_cols].values).all(axis=1)
            
            if duplicate_mask.any():
                existing_df.loc[duplicate_mask, "lisi_score"] = val_lisi
                existing_df.to_csv(lisi_log, mode="w", header=True, index=False)
                print(f"Updated existing entry in {lisi_log}")
            else:
                new_row.to_csv(lisi_log, mode="a", header=False, index=False)
                print(f"Appended new entry to {lisi_log}")
        else:
            new_row.to_csv(lisi_log, mode="w", header=True, index=False)
            print(f"Created new log at {lisi_log}")

    # ── Time-supervised training ───────────────────────────────────────────
    if predict == 'time':
        latent_means = np.load(embeddings_path)
        adata.obsm["X_cvae"] = latent_means

        (data_train, data_train_label,
         data_val_t,  data_val_label,
         _,   _) = split_time_trainval(
            adata, train_species, seed, nsubsample=10000)

        print(f'time pred training set shape: {data_train.shape}')
        print(f'time pred validation set shape: {data_val_t.shape}')
        batch_size_time = 512
        batch_size_val  = 512

        # shared folder for all temporal models derived from this CVAE embedding
        time_model_outdir = Path(f'{outdir}/time_prediction_models/{model_name}')
        time_model_outdir.mkdir(parents=True, exist_ok=True)

        # JSON record of best hyperparameters — keyed by train_species
        best_model_record = time_model_outdir / f'best_model_{train_species}.json'

        if best_model_record.exists():
            # ── reload best hparams from previous sweep ────────────────────
            print('Best temporal model record found. Skipping sweep...')
            with open(best_model_record) as f:
                cfg = json.load(f)
            lr_t = cfg['lr']
            nl   = cfg['nlayer']
            ed   = cfg['embed_dim']
            print(f'  Best hparams: lr={lr_t}  nlayer={nl}  embed_dim={ed}')

        else:
            # ── hyperparameter sweep ───────────────────────────────────────
            learning_rate_list = [0.1, 0.01, 0.001, 0.0001]
            nlayer_list        = [2, 3, 4]
            embed_dim_list     = [50, 100, 200, 400, 800]

            mse_list = []
            dic      = {}
            i        = 0
            
            for lr_t in learning_rate_list:
                for nl in nlayer_list:
                    for ed in embed_dim_list:
                        dic[i] = [lr_t, nl, ed]
                        print(f'Sweep [{i}]: lr={lr_t}  nlayer={nl}  embed_dim={ed}')

                        sweep_dir = time_model_outdir / f'sweep_{lr_t}_{nl}_{ed}'
                        predictor = TimePredictor(
                            input_dim=data_train.shape[1], embed_dim=ed,
                            nlayer=nl, output_model=str(sweep_dir),
                            learning_rate_x=lr_t, device=device)
                        predictor.train(
                            data_train, data_val_t,
                            data_train_label, data_val_label,
                            output_model=str(sweep_dir),
                            batch_size=batch_size_time, my_epochs=200, patience=20)

                        nbatch_val = data_val_t.shape[0] // batch_size_val
                        mse_list_i = [
                            predictor.get_losses_rna(
                                data_val_t[batch_size_val * b:
                                           min(batch_size_val * (b + 1), data_val_t.shape[0])],
                                data_val_label[batch_size_val * b:
                                               min(batch_size_val * (b + 1), data_val_t.shape[0])])
                            for b in range(nbatch_val)]
                        mse = sum(mse_list_i) / len(mse_list_i)
                        print(f'  Val MSE: {mse:.6f}')
                        mse_list.append(mse)
                        i += 1

            # pick best and save record
            best_i       = int(np.argmin(mse_list))
            lr_t, nl, ed = dic[best_i]
            print(f'Best hparams: lr={lr_t}  nlayer={nl}  embed_dim={ed}  '
                  f'(val MSE={mse_list[best_i]:.6f})')
            cfg = {
                'lr':        lr_t,
                'nlayer':    nl,
                'embed_dim': ed,
                'input_dim': int(data_train.shape[1]),
            }
            with open(best_model_record, 'w') as f:
                json.dump(cfg, f, indent=2)
            print(f'  Saved best model record -> {best_model_record}')

        # ── train / reload best model into fixed path ──────────────────────
        best_model_dir = time_model_outdir / 'best_model'
        print(f'Training/loading best temporal model at {best_model_dir}...')
        predictor = TimePredictor(
            input_dim=data_train.shape[1], embed_dim=ed,
            nlayer=nl, output_model=str(best_model_dir),
            learning_rate_x=lr_t, device=device)
        predictor.train(
            data_train, data_val_t,
            data_train_label, data_val_label,
            output_model=str(best_model_dir),
            batch_size=batch_size_time, my_epochs=200, patience=20)

        # ensure input_dim is always up to date in the record
        cfg['input_dim'] = int(data_train.shape[1])
        with open(best_model_record, 'w') as f:
            json.dump(cfg, f, indent=2)

        # ── predict and save ───────────────────────────────────────────────
        adata_target = adata[adata.obs.species == target_species]
        predictions  = predictor.predict_time(adata_target.obsm["X_cvae"])
        adata_target.obs['pred_time'] = predictions

        out_csv = f'{outdir}/{model_name}_pred_time_{target_species}.txt'
        adata_target.obs.to_csv(out_csv, index=False, sep='\t', float_format="%.4f")
        print(f'  Saved predictions -> {out_csv}')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Optional app description')

    print("Command being run:")
    print("python", " ".join(sys.argv))

    parser.add_argument('--input_h5ad',            type=str,   help='input_h5ad')
    parser.add_argument('--learning_rate',         type=float, help='learning_rate',                           default=0.001)
    parser.add_argument('--predict',               type=str,   help='options: train or time',                  default='')
    parser.add_argument('--nlayer',                type=int,   help='nlayer',                                  default=3)
    parser.add_argument('--batch_size',            type=int,   help='batch size',                              default=128)
    parser.add_argument('--my_epochs',             type=int,   help='maximum number of epochs',                default=500)
    parser.add_argument('--dropout_rate',          type=float, help='dropout_rate',                            default=0.1)
    parser.add_argument('--embed_dim',             type=int,   help='embed_dim',                               default=25)
    parser.add_argument('--patience',              type=int,   help='patience',                                default=25)
    parser.add_argument('--hidden_size',           type=int,   help='hidden_size',                             default=512)
    parser.add_argument('--target_species',        type=str,   help='target species',                          default='')
    parser.add_argument('--train_species',         type=str,   help='species to train',                        default='')
    parser.add_argument('--target_time',           type=str,   help='target time',                             default='')
    parser.add_argument('--group',                 type=str,   help='obs column for group info',               default='')
    parser.add_argument('--batch_colname',         type=str,   help='column name of obs that correspond to batch information)', default='batch')
    parser.add_argument('--dis',                   type=str,   help='"dis" to use discriminator, "" otherwise', default='')
    parser.add_argument('--discriminator_weight',  type=float, help='weight of discriminator loss in VAE generator step', default=1.0)
    parser.add_argument('--time_label',            type=str,   help='obs column for time label',                default='time')
    parser.add_argument('--celltype_col',          type=str, help='obs column for cell type',                   default='cell_type')
    parser.add_argument('--seed',                  type=int, help='random seed',                                default='101')

    args = parser.parse_args()

    main(args)
