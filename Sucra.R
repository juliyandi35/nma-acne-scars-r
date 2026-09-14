library(gemtc)
library(igraph)
library(rjags)

# Gunakan data kamu langsung
data_ab <- blocker$data.ab

network <- mtc.network(data.ab = data_ab)

plot(network, layout = layout.fruchterman.reingold)

model <- mtc.model(
  network,
  linearModel = "fixed",          # bisa diganti "random" kalau mau
  likelihood = "binom",
  link = "logit",
  n.chain = 4
)

mcmc <- mtc.run(
  model,
  n.adapt = 5000,
  n.iter = 100000,
  thin = 10
)

rp <- rank.probability(mcmc)

sucra_res <- sucra(rp)

sucra_res

plot(sucra_res)
