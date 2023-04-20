library(io)

# two studies were both done on GPL10558
res.g1 <- qread("../GSE49135/GSE49135_limma.rds");
res.g2 <- qread("../GSE62061/GSE62061_limma.rds");

hist(res.g1[[1]]$t, breaks=100)
hist(res.g2[[1]]$t, breaks=100)

filter_rank <- function(res, K) {
	r <- rank(res$t);
	idx <- r <= K | r > nrow(res) - K;
	res[idx, ]
}

# filter for top k genes
K <- 3000;
res.g1.f <- lapply(res.g1, filter_rank, K=K);
res.g2.f <- lapply(res.g2, filter_rank, K=K);

# match probes so that they are in the same order across studies
# NB  if comparing across different platforms,
#     we should match based on ensembl IDs
probes <- intersect(
	Reduce(function(x, y) intersect(x, y), lapply(res.g1.f, rownames)),
	Reduce(function(x, y) intersect(x, y), lapply(res.g2.f, rownames))
);
stopifnot(length(probes) > 0);
res.g1.m <- lapply(res.g1, function(res) res[probes, ]);
res.g2.m <- lapply(res.g2, function(res) res[probes, ]);

# check that the probes are in the same order
for (i in 1:length(res.g1.m)) {
	for (j in 1:length(res.g2.m)) {
		stopifnot(rownames(res.g1.m[[i]]) == rownames(res.g2.m[[j]]))
	}
}

# ---

smooth_scatter <- function(x, y, cut=K, ...) {
	xlim <- range(x);
	ylim <- range(y);
	lim <- c(min(xlim[1], ylim[1]), max(xlim[2], ylim[2]));
	if (length(x) > 1e3) {
		smoothScatter(x, y, xlim=lim, ylim=lim, nbin=512, ...);
		points(x, y, pch='.');
	} else {
		plot(x, y, xlim=lim, ylim=lim, ...);
	}
	abline(a=0, b=1, col="grey30");
	# abline(
	# 	h = c(cut, length(y) - cut),
	# 	v = c(cut, length(x) - cut),
	# 	col = "grey30"
	# );
	cor(x,  y)
}

# ideal correlation for the same cell line (i.e. with itself)
smooth_scatter(res.g1.m[[1]]$t, res.g1.m[[1]]$t);

# cell line 1 in study 1 vs. cell line 1 in study 2
smooth_scatter(res.g1.m[[1]]$t, res.g2.m[[1]]$t);

# cell line 1 vs. 2 in study 2
smooth_scatter(res.g2.m[[1]]$t, res.g2.m[[2]]$t);

# cell line 1 vs. 3 in study 2
smooth_scatter(res.g2.m[[1]]$t, res.g2.m[[3]]$t);

# cell line 1 vs. 4 in study 2
smooth_scatter(res.g2.m[[1]]$t, res.g2.m[[4]]$t);

# ---

# add rank of t statistic to results
add_rank <- function(res) {
	res$rank.t <- rank(res$t);
	res
}

res.g1.m <- lapply(res.g1.m, add_rank);
res.g2.m <- lapply(res.g2.m, add_rank);

# between-study comparisons

# cell line 1 in study 1 vs. cell line 1 in study 2
smooth_scatter(res.g1.m[[1]]$rank.t, res.g2.m[[1]]$rank.t);

# cell line 1 in study 1 vs. cell line 2 in study 2
smooth_scatter(res.g1.m[[1]]$rank.t, res.g2.m[[2]]$rank.t);

# cell line 1 in study 1 vs. cell line 2 in study 2
smooth_scatter(res.g1.m[[1]]$rank.t, res.g2.m[[3]]$rank.t);

# cell line 1 in study 1 vs. cell line 2 in study 2
smooth_scatter(res.g1.m[[1]]$rank.t, res.g2.m[[4]]$rank.t);

# within-study comparisons

# cell line 1 vs. 2 in study 2
smooth_scatter(res.g2.m[[1]]$rank.t, res.g2.m[[2]]$rank.t);

# cell line 1 vs. 3 in study 2
smooth_scatter(res.g2.m[[1]]$rank.t, res.g2.m[[3]]$rank.t);

# cell line 1 vs. 4 in study 2
smooth_scatter(res.g2.m[[1]]$rank.t, res.g2.m[[4]]$rank.t);

