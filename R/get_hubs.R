get_cor_cut <- function(data_all, MR_E, cut_Data, RR = 10) {
  ff <- c()
  for(i in 1:RR) {
    sample_data <- data_all[sample(1:nrow(data_all), cut_Data), ]
    cor_C0 <- cor(t(sample_data))
    cor_C0[is.na(cor_C0)] <- 0
    MRHCA_output <- MRHCA(cor_C0, MR_E, step_size0 = 50, p_sig = 0.05, Hit_score_cutoff = 100)
    tg_ids<-which(MRHCA_output [[1]][,1]==1)
    cor_all_c <- c()
    for(j in 1:length(tg_ids)) {
      x <- cor(data_all[tg_ids[j], ], t(data_all))
      x[is.na(x)] <- 0
      cor_all_c <- c(cor_all_c, sort(x, decreasing = T)[MRHCA_output[[1]][tg_ids[j], 2]])
    }
    ff <- c(ff, mean(cor_all_c))
  }
  cor_cut<-median(ff)
}

GetCorCut <- function(data_all, cut_Data = 1000, RR = 10) {
  sample_data_r <- data_all[sample(1:nrow(data_all), cut_Data), ]
  for(i in 1:nrow(sample_data_r)) {
    sample_data_r[i, ] <- sample(sample_data_r[i, ])
  }
  
  MR_E <- Empirical_null_distribution_MR(sample_data_r, Rounds = 1000)
  cor_cut <- get_cor_cut(data_all, MR_E, cut_Data, RR)
}

GetHubs <- function(x, tn.p = 1, k = 500, RR = 10) {
  stopifnot(tn.p > 0)
  stopifnot(tn.p <= 1)
  stopifnot(k <= nrow(x))
  
  data_all <- x[which(apply(x != 0, 1, sum) > 5), ]
  
  mrs <- GetHubInfos(data_all, tn.p, k)
  
  hubs <- matrix(unlist(mrs[[1]], use.names = FALSE), ncol = 3, byrow = TRUE)
  mr_id <- mrs[[2]] + 1
  mr.em <- mrs[[3]]
  hub.indices <- which(hubs[, 1] > 0)

  if(nrow(data_all) > k) {
    cor_cut <- GetCorCut(data_all, cut_Data = k, RR = RR)
    
    tg_large_hubs <- which((hubs[, 1] == 1) & (hubs[, 2] == 0))
    for(i in 1:length(tg_large_hubs)) {
      cor_all_c <- cor(data_all[tg_large_hubs[i], ], t(data_all))[1, ]
      hubs[tg_large_hubs[i], 3] <- sum(cor_all_c > cor_cut)
    }
  }
  
  selected <- sapply(hub.indices, function(index){
    count <- hubs[index, 2] + hubs[index, 3]
    mr <- mr.em[[index]]
    mr < sort(mr, partial = count)[count]
  })
  
  if (length(selected) > 0) rownames(selected) <- rownames(data_all)
  t(selected)
}

FixHubsData <- function(data, id, em, cor_cut, k = 500, RR = 10, filter = 20) {
  tg_large_hubs <- which((id[, 1] == 1) & (id[, 2] == 0))
  for(i in 1:length(tg_large_hubs)) {
    cor_all_c <- cor(data[id[tg_large_hubs[i], 4] + 1, ], t(data))[1, ]
    cor_all_c[is.na(cor_all_c)] <- 0
    id[tg_large_hubs[i], 3] <- sum(cor_all_c > cor_cut)
  }
  
  id_selected <- id[which(id[, 2] + id[, 3] > filter), ]
  em_selected <- em[which(id[, 2] + id[, 3] > filter), ]
  
  selected <- sapply(1:dim(id_selected)[1], function(index){
    count <- id_selected[index, 2] + id_selected[index, 3]
    mr <- unlist(em_selected[index, ])
    mr < sort(mr, partial = count)[count]
  })
  
  if (length(selected) > 0) rownames(selected) <- rownames(data)
  t(selected)
}

FixHubs <- function(datafile, file, emfile, k = 500, RR = 10, filter = 20) {
  col.names <- unlist(strsplit(scan(datafile, what = 'character', nlines = 1, sep = "\n", quiet = TRUE), split = '\t'))[-1]
  data <- bigmemory::as.matrix(bigmemory::read.big.matrix(file = datafile, sep = "\t", skip = 1, header = FALSE, col.names = col.names, has.row.names = TRUE, type = "double"))
  id <- data.table::fread(file = file, sep = "\t")
  em <- data.frame(data.table::fread(file = emfile, sep = "\t"), row.names = 1)
  stopifnot(dim(id)[1] == dim(em)[1])
  stopifnot(all(id[, 4] == row.names(em)))
  cor_cut <- GetCorCut(data, cut_Data = k, RR = RR)
  FixHubsData(data, id, em, cor_cut, k, RR, filter)
}

FixHubsWithCorCut <- function(datafile, file, emfile, cor_cut, k = 500, RR = 10, filter = 20) {
  col.names <- unlist(strsplit(scan(datafile, what = 'character', nlines = 1, sep = "\n", quiet = TRUE), split = '\t'))[-1]
  data <- bigmemory::as.matrix(bigmemory::read.big.matrix(file = datafile, sep = "\t", skip = 1, header = FALSE, col.names = col.names, has.row.names = TRUE, type = "double"))
  print("data loaded")
  id <- data.table::fread(file = file, sep = "\t")
  print("id loaded")
  em <- data.frame(data.table::fread(file = emfile, sep = "\t"), row.names = 1)
  print("em loaded")
  stopifnot(dim(id)[1] == dim(em)[1])
  stopifnot(all(id[, 4] == row.names(em)))
  if(missing(cor_cut)) cor_cut <- GetCorCut(data, cut_Data = k, RR = RR)
  print(cor_cut)
  FixHubsData(data, id, em, cor_cut, k, RR, filter)
}