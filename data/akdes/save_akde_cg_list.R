###
akde_path <- "./data/akdes/muddyfoot_roach_akdes/akde_cg/"


# List all RDS files in the folder
rds_files <- list.files(path = akde_path, pattern = "\\.rds$", full.names = TRUE)

# Read all RDS files into a list
rds_list <- lapply(rds_files, readRDS)

# Read all RDS files into a list
names(rds_list) <- basename(rds_files)
names(rds_list) <- sub("\\_akde_cg.rds$", "", names(rds_list))

# Print the names of the loaded RDS files
print(names(rds_list))

roach_akdes_cg_list <- rds_list

summary(roach_akdes_cg_list$F59704)


#save
saveRDS(roach_akdes_cg_list, paste0(akde_path, "roach_akdes_cg_list.rds"))
