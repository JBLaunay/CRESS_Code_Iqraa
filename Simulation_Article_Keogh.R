### Chargement des bibliothèques
install.packages("survival")
install.packages("dplyr")
install.packages("ggplot2")


library(survival)
library(dplyr)
library(ggplot2)


# Test pour vérifier que Git fonctionne bien

########################## DEBUT DE LA SIMULATION #####################################

### Paramètres de la simulation
set.seed(123)  # Pour la reproductibilité
n <- 5000       # Nombre d'individus
unite_de_temps <- 4       # Nombre de périodes
nombre_de_mesures <- 10      # Nombre de mesures par unité de temps
k_max <- unite_de_temps * nombre_de_mesures


# Paramètres pour la simulation
delta_0 <- 0
delta_L <- 0.8
delta_A <- -1
delta_T <- 0.1
gamma_0 <- -1
gamma_L <- 0.5
gamma_A <- 1e5
alpha_0 <- 0.2
alpha_A <- -0.04
alpha_L <- 0.015
alpha_U <- 0.015


# Paramètres pour la simulation modifiés pour tenir compte de la dilatation temporelle

delta_0_corrigé <- delta_0/nombre_de_mesures
delta_L_corrigé<- delta_L^(1/nombre_de_mesures)
delta_A_corrigé <- delta_A/nombre_de_mesures
delta_T_corrigé <- (1-delta_L^(1/nombre_de_mesures))/((1-delta_L) * nombre_de_mesures)
delta_U_corrigé <- 1/nombre_de_mesures
standard_deviation_corrigé <- 1/nombre_de_mesures^(1/2)

### Simulation des données
print("Initialisation des matrices L et A...")
U <- rnorm(n, mean = 0, sd = sqrt(0.1))  # Variable non observée
L_0 <- rnorm(n, mean = U, sd = 1)
L <- matrix(NA, nrow = n, ncol = k_max + 1)
A <- matrix(0, nrow = n, ncol = k_max + 1)
T <- rep(NA, n)

L[, 1] <- L_0
A[, 1] <- ifelse(runif(n) < plogis(gamma_0 + gamma_L * L[, 1]), 1, 0)  # Initialisation de A pour k = 1

for (k in 1:k_max) {
  L[, k + 1] <- rnorm(n, mean = delta_0_corrigé + delta_L_corrigé * L[, k] + delta_A_corrigé * A[, k] + delta_T_corrigé * k + delta_U_corrigé * U, sd = standard_deviation_corrigé)
  logit_A <- gamma_0 + gamma_L * L[, k + 1] + gamma_A * A[, k]
  A[, k + 1] <- ifelse(runif(n)^(1/nombre_de_mesures) < plogis(logit_A), 1, 0)
}


print("Matrices L et A initialisées.")
print(head(L))  # Affiche les premières lignes de L
print(head(A))  # Affiche les premières lignes de A


### Simulation des temps de survie
hazard <- function(k, A, L, U) {
  alpha_0 + alpha_A * A + alpha_L * L + alpha_U * U
}

# Création du tableau des risques instantanés
risques_instantanés <- matrix(NA, nrow = n, ncol = k_max + 1)
for (k in 1:(k_max + 1)) {
  risques_instantanés[, k] <- hazard(k, A[, k], L[, k], U)
}

# Création du tableau des risques cumulés
risques_cumulés <- matrix(NA, nrow = n, ncol = k_max + 1)
for (j in 1:(k_max + 1)) {
  risques_cumulés[, j] <- rowSums(risques_instantanés[, 1:j, drop = FALSE])
}

# Création du tableau des probabilités de survie
proba_de_survie <- exp(-risques_cumulés/nombre_de_mesures)

# Affichage des résultats
print(head(risques_instantanés))  # Affiche les premières lignes des risques instantanés
print(head(risques_cumulés))      # Affiche les premières lignes des risques cumulés
print("Proba de survie simulées.")
print(head(proba_de_survie))      # Affiche les premières lignes des probabilités de survie

# Tirer n nombres uniformément entre 0 et 1
nombres_uniformes <- runif(n)

# # Création du tableau des temps de survie sans censure
# temps_de_survie <- matrix(NA, nrow = n, ncol = k_max + 1)
# for (j in 1:(k_max + 1)) {
#   temps_de_survie[, j] <- ifelse(proba_de_survie[, j] < nombres_uniformes, 1, 0)
# }

# Création du tableau des temps de survie avec censure
temps_de_survie <- matrix(NA, nrow = n, ncol = k_max + 1)
for (j in 1:(k_max + 1)) {
  temps_de_survie[, j] <- ifelse(proba_de_survie[, j] < nombres_uniformes, 1,
                                 ifelse(runif(n) < 0.1, -1, 0))  # 10% de censure
}

# Identifier les individus traités et non traités
traitement <- A[, k_max + 1]


### Tracé des courbes de survie

# # Calculer la proportion de patients encore vivants pour chaque période et chaque groupe sans censure
# proportion_vivants <- data.frame(
#   temps = rep(1:(k_max + 1), 2),
#   groupe = rep(c("Traités", "Non traités"), each = k_max + 1),
#   proportion = c(
#     sapply(1:(k_max + 1), function(j) mean(temps_de_survie[traitement == 1, j] == 0)),
#     sapply(1:(k_max + 1), function(j) mean(temps_de_survie[traitement == 0, j] == 0))
#   )
# )

# Création des données pour l'estimateur de Kaplan-Meier
survie_data <- data.frame(
  time = apply(temps_de_survie, 1, function(row) min(which(row != 0))),
  status = apply(temps_de_survie, 1, function(row) ifelse(row[min(which(row != 0))] == 1, 1, 0)),
  groupe = ifelse(traitement == 1, "Traités", "Non traités")
)


# Estimation des courbes de survie de Kaplan-Meier
print(km_fit$time)
print(km_fit$surv)
print(km_fit$strata)
km_fit <- survfit(Surv(time, status) ~ groupe, data = survie_data)

# Calcul de L en moyenne
km_data <- data.frame(
  time = rep(km_fit$time, each = length(km_fit$strata)),
  surv = rep(km_fit$surv, each = length(km_fit$strata)),
  strata = factor(rep(km_fit$strata, each = length(km_fit$time)), labels
                  = c("Non traités", "Traités"))
)

ggplot(km_data, aes(x = time, y = surv, color = strata)) +
  geom_step() +
  labs(title = "Courbes de survie de Kaplan-Meier",
       x = "Temps",
       y = "Probabilité de survie") +
  theme_minimal()

km_fit_traités <- survfit(Surv(time, status) ~ 1, data = survie_data[survie_data$groupe=="Traités",])
km_fit_non_traités <- survfit(Surv(time, status) ~ 1, data = survie_data[survie_data$groupe=="Non traités",])


km_data_traités <- data.frame(
  time = km_fit_traités$time,
  surv = km_fit_traités$surv,
  strata = "Traités"
)

km_data_non_traités <- data.frame(
  time = km_fit_non_traités$time,
  surv = km_fit_non_traités$surv,
  strata = "Non traités"
)

ggplot(km_data_traités, aes(x = time, y = surv)) +
  geom_step() +
  labs(title = "Courbe de survie de Kaplan-Meier pour les traités",
       x = "Temps",
       y = "Probabilité de survie") +
  theme_minimal()

ggplot(km_data_non_traités, aes(x = time, y = surv)) +
  geom_step() +
  labs(title = "Courbe de survie de Kaplan-Meier pour les non traités",
       x = "Temps",
       y = "Probabilité de survie") +
  theme_minimal()
