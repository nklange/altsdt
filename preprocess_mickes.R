# General information ----------------

# preprocessing for SB2020 and SB2022 straightforward based on authors' osf code
# for MWW2007, preprocessing below as it includes some binning of confidence
# ratings for different versions of the data


library(dplyr)
library(tidyr)
library(magrittr)
library(purrr)
library(readr)
library(tibble)
library(readxl)

Exp1<-read_excel("EmpiricalData/MWWData/Data for David Kellen.xls", sheet = "Exp1") %>% mutate(exp = "MWW2007_e1") %>% mutate(id = paste0(exp,"_",id))
Exp2<-read_excel("EmpiricalData/MWWData/Data for David Kellen.xls", sheet = "Exp2") %>% mutate(exp = "MWW2007_e2") %>% mutate(id = paste0(exp,"_",id))
Exp3<-read_excel("EmpiricalData/MWWData/Data for David Kellen.xls", sheet = "Exp3") %>% mutate(exp = "MWW2007_e3") %>% mutate(id = paste0(exp,"_",id))





Exp2n_twentypoint <- Exp2 %>% mutate(response = case_when(response %in% c(1:5) ~ 1,
                                              response %in% c(6:10) ~ 2,
                                              response %in% c(11:15) ~ 3,
                                              response %in% c(16:20) ~ 4,
                                              response %in% c(21:25) ~ 5,
                                              response %in% c(26:30) ~ 6,
                                              response %in% c(31:35) ~ 7,
                                              response %in% c(36:40) ~ 8,
                                              response %in% c(41:45) ~ 9,
                                              response %in% c(46:50) ~ 10,
                                              response %in% c(51:55) ~ 11,
                                              response %in% c(56:60) ~ 12,
                                              response %in% c(61:65) ~ 13,
                                              response %in% c(66:70) ~ 14,
                                              response %in% c(71:75) ~ 15,
                                              response %in% c(76:80) ~ 16,
                                              response %in% c(81:85) ~ 17,
                                              response %in% c(86:90) ~ 18,
                                              response %in% c(91:95) ~ 19,
                                              response %in% c(96:99) ~ 20
                                              ))

Exp2n_sixpoint <- Exp2 %>% mutate(response = case_when(response %in% c(1:16) ~ 1,
                                                          response %in% c(17:32) ~ 2,
                                                          response %in% c(33:49) ~ 3,
                                                          response %in% c(50:66) ~ 4,
                                                          response %in% c(67:82) ~ 5,
                                                          response %in% c(83:99) ~ 6
))

Exp1n_sixpoint <- Exp1 %>% mutate(response = case_when(response %in% c(1:4) ~ 1,
                                                       response %in% c(5:7) ~ 2,
                                                       response %in% c(8:10) ~ 3,
                                                       response %in% c(11:13) ~ 4,
                                                       response %in% c(14:16) ~ 5,
                                                       response %in% c(17:20) ~ 6
))





MWW2007 <- bind_rows(Exp1n_sixpoint,Exp2n_sixpoint,Exp3) %>%
  mutate(oldnew = ifelse(oldnew=="y","Old","New")) %>%
  filter(oldnew %in% c("Old","New"))

saveRDS(MWW2007,"EmpiricalData/MWWData/MWW2007_sixpoint_preprocesseddata.rds")

