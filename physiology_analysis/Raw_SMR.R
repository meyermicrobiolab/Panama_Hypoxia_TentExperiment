library(ggplot2)
library(respR)
library(viridis)
library(ggplot2)
library(respR)
library(ggpubr)
library(respirometry)
library(vegan)
library(cowplot)
library(mvabund) 
library(viridis)
library(tidyverse)
library(car) # homogenity of variances
library(lattice)  
library(mgcv) #Load the mgcv package

remove(list = ls())

# top ---------------------------------------------------------------------
df <- read.csv(file = "Raw_SMR_tent.csv")
df <- na.omit(df)
str(df)

# respiration objects for Aten
df_aten <- df[which(df$spp == 'aten'),]
df_sid <- df[which(df$spp != 'aten'),]

# agaricia
aten24 <- df_aten[which(df_aten$ID == '24'),]
aten2 <- df_aten[which(df_aten$ID == '2'),]
aten3 <- df_aten[which(df_aten$ID == '3'),]
aten11 <- df_aten[which(df_aten$ID == '11'),]
aten21 <- df_aten[which(df_aten$ID == '21'),]
aten13 <- df_aten[which(df_aten$ID == '13'),]
aten5 <- df_aten[which(df_aten$ID == '5'),]
aten8 <- df_aten[which(df_aten$ID == '8'),]
aten22 <- df_aten[which(df_aten$ID == '22'),]
aten_blank <- df_aten[which(df_aten$ID == 'blank'),]
aten17 <- df_aten[which(df_aten$ID == '17'),]
aten26 <- df_aten[which(df_aten$ID == '26'),]
aten10 <- df_aten[which(df_aten$ID == '10'),]
aten28 <- df_aten[which(df_aten$ID == '28'),]

# siderastrea
sid_blank <- df_sid[which(df_sid$ID == 'blank1'),]

sid7 <- df_sid[which(df_sid$ID == '7'),]
sid15 <- df_sid[which(df_sid$ID == '15'),]
sid16 <- df_sid[which(df_sid$ID == '16'),]
sid9 <- df_sid[which(df_sid$ID == '9'),]
sid6 <- df_sid[which(df_sid$ID == '6'),]
sid12 <- df_sid[which(df_sid$ID == '12'),]
sid4 <- df_sid[which(df_sid$ID == '4'),]
sid31 <- df_sid[which(df_sid$ID == '31'),]
sid32 <- df_sid[which(df_sid$ID == '32'),]
sid25 <- df_sid[which(df_sid$ID == '25'),]
sid14 <- df_sid[which(df_sid$ID == '14'),]
sid23 <- df_sid[which(df_sid$ID == '23'),]
sid29 <- df_sid[which(df_sid$ID == '29'),]
sid27 <- df_sid[which(df_sid$ID == '27'),]
sid19 <- df_sid[which(df_sid$ID == '19'),]
sid20 <- df_sid[which(df_sid$ID == '20'),]


bg1 <- calc_rate.bg(aten_blank, time = 9, oxygen = 11, plot = TRUE)
bg1 # unit-less rates returned, yippee!

# 24 ---------------------------------------------------------------------

inspect(aten24, time = 9, oxygen = 11, plot = TRUE)
aten24r <- data.frame(aten24$delta_t, aten24$Value)

# aten24 <- subset_data(aten24, from = 205.688, to = 163.914, by = "O2") # subsetting 

rate_aten24 <- calc_rate(aten24r, plot=TRUE)
rate_aten24_adj <- adjust_rate(rate_aten24r, bg1)

R_aten24 <- convert_rate(rate_aten24_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                           volume = 0.1379, 
                           area = .023669576, 
                           S = aten24$Salinity, 
                           t = aten24$Temp)
R_aten24

R_aten24_abs <- convert_rate(rate_aten24_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                         volume = 0.1379, 
                         #area = .023669576, 
                         S = aten24$Salinity, 
                         t = aten24$Temp)
R_aten24_abs

# 2 ---------------------------------------------------------------------
inspect(aten2, time = 9, oxygen = 11, plot = TRUE)
aten2r <- data.frame(aten2$delta_t, aten2$Value)
rate_aten2 <- calc_rate(aten2r, plot=TRUE)
rate_aten2_adj <- adjust_rate(rate_aten2, bg1)

R_aten2 <- convert_rate(rate_aten2_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                         volume =0.1361, 
                         area = 0.02814924, 
                         S = aten2$Salinity, t = aten2$Temp)
R_aten2

R_aten2_abs <- convert_rate(rate_aten2_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                             volume = 0.1361, 
                             S = aten2$Salinity, t = aten2$Temp)
R_aten2_abs

# 3 ---------------------------------------------------------------------
inspect(aten3, time = 9, oxygen = 11, plot = TRUE)
aten3r <- data.frame(aten3$delta_t, aten3$Value)
rate_aten3 <- calc_rate(aten3r, plot=TRUE)
rate_aten3_adj <- adjust_rate(rate_aten3, bg1)

R_aten3 <- convert_rate(rate_aten3_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume =0.1412, 
                        area = 0.020688178, 
                        S = aten3$Salinity, t = aten3$Temp)
R_aten3

R_aten3_abs <- convert_rate(rate_aten3_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1412, 
                            S = aten3$Salinity, t = aten3$Temp)
R_aten3_abs

# 11 ---------------------------------------------------------------------
inspect(aten11, time = 9, oxygen = 11, plot = TRUE)
aten11r <- data.frame(aten11$delta_t, aten11$Value)
rate_aten11 <- calc_rate(aten11r, plot=TRUE)
rate_aten11_adj <- adjust_rate(rate_aten11, bg1)

R_aten11 <- convert_rate(rate_aten11_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume =0.1392, 
                        area = 0.016435524, 
                        S = aten11$Salinity, t = aten11$Temp)
R_aten11

R_aten11_abs <- convert_rate(rate_aten11_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1392, 
                            S = aten11$Salinity, t = aten11$Temp)
R_aten11_abs

# 21 ---------------------------------------------------------------------
inspect(aten21, time = 9, oxygen = 11, plot = TRUE)
aten21r <- data.frame(aten21$delta_t, aten21$Value)
rate_aten21 <- calc_rate(aten21r, plot=TRUE)
rate_aten21_adj <- adjust_rate(rate_aten21, bg1)

R_aten21 <- convert_rate(rate_aten21_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume =0.1352, 
                        area = 0.062140204, 
                        S = aten21$Salinity, t = aten21$Temp)
R_aten21

R_aten21_abs <- convert_rate(rate_aten21_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1352, 
                            S = aten21$Salinity, t = aten21$Temp)
R_aten21_abs

# 13 ---------------------------------------------------------------------
inspect(aten13, time = 9, oxygen = 11, plot = TRUE)
aten13r <- data.frame(aten13$delta_t, aten13$Value)
rate_aten13 <- calc_rate(aten13r, plot=TRUE)
rate_aten13_adj <- adjust_rate(rate_aten13, bg1)

R_aten13 <- convert_rate(rate_aten13_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume =0.1338, 
                        area = 0.062140204, 
                        S = aten13$Salinity, t = aten13$Temp)
R_aten13

R_aten13_abs <- convert_rate(rate_aten13_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1338, 
                            S = aten13$Salinity, t = aten13$Temp)
R_aten13_abs

# 5 ---------------------------------------------------------------------
inspect(aten5, time = 9, oxygen = 11, plot = TRUE)
aten5r <- data.frame(aten5$delta_t, aten5$Value)
rate_aten5 <- calc_rate(aten5r, plot=TRUE)
rate_aten5_adj <- adjust_rate(rate_aten5, bg1)

R_aten5 <- convert_rate(rate_aten5_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume =0.1379, 
                        area = 0.043328642, 
                        S = aten5$Salinity, t = aten5$Temp)
R_aten5

R_aten5_abs <- convert_rate(rate_aten5_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1379, 
                            S = aten5$Salinity, t = aten5$Temp)
R_aten5_abs

# 8 ---------------------------------------------------------------------
inspect(aten8, time = 9, oxygen = 11, plot = TRUE)
aten8r <- data.frame(aten8$delta_t, aten8$Value)
rate_aten8 <- calc_rate(aten8r, plot=TRUE)
rate_aten8_adj <- adjust_rate(rate_aten8, bg1)

R_aten8 <- convert_rate(rate_aten8_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume =0.1367, 
                        area = 0.04805045, 
                        S = aten8$Salinity, t = aten8$Temp)
R_aten8

R_aten8_abs <- convert_rate(rate_aten8_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1367, 
                            S = aten8$Salinity, t = aten8$Temp)
R_aten8_abs

# 22 ---------------------------------------------------------------------
inspect(aten22, time = 9, oxygen = 11, plot = TRUE)
aten22r <- data.frame(aten22$delta_t, aten22$Value)
rate_aten22 <- calc_rate(aten22r, plot=TRUE)
rate_aten22_adj <- adjust_rate(rate_aten22, bg1)

R_aten22 <- convert_rate(rate_aten22_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume = 0.1385, 
                        area = 0.034550922, 
                        S = aten22$Salinity, t = aten22$Temp)
R_aten22

R_aten22_abs <- convert_rate(rate_aten22_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1385, 
                            S = aten22$Salinity, t = aten22$Temp)
R_aten22_abs

# 17 ---------------------------------------------------------------------
inspect(aten17, time = 9, oxygen = 11, plot = TRUE)
aten17r <- data.frame(aten17$delta_t, aten17$Value)
rate_aten17 <- calc_rate(aten17r, plot=TRUE)
rate_aten17_adj <- adjust_rate(rate_aten17, bg1)

R_aten17 <- convert_rate(rate_aten17_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume = 0.1387, 
                        area = 0.02686285, 
                        S = aten17$Salinity, t = aten17$Temp)
R_aten17

R_aten17_abs <- convert_rate(rate_aten17_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1387, 
                            S = aten17$Salinity, t = aten17$Temp)
R_aten17_abs

# 26 ---------------------------------------------------------------------
inspect(aten26, time = 9, oxygen = 11, plot = TRUE)
aten26r <- data.frame(aten26$delta_t, aten26$Value)
rate_aten26 <- calc_rate(aten26r, plot=TRUE)
rate_aten26_adj <- adjust_rate(rate_aten26, bg1)

R_aten26 <- convert_rate(rate_aten26_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume = 0.1360, 
                        area = 0.02035523, 
                        S = aten26$Salinity, t = aten26$Temp)
R_aten26

R_aten26_abs <- convert_rate(rate_aten26_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1360, 
                            S = aten26$Salinity, t = aten26$Temp)
R_aten26_abs

# 10 ---------------------------------------------------------------------
inspect(aten10, time = 9, oxygen = 11, plot = TRUE)
aten10r <- data.frame(aten10$delta_t, aten10$Value)
rate_aten10 <- calc_rate(aten10r, plot=TRUE)
rate_aten10_adj <- adjust_rate(rate_aten10, bg1)

R_aten10 <- convert_rate(rate_aten10_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume = 0.1357, 
                        area = 0.052015558, 
                        S = aten10$Salinity, t = aten10$Temp)
R_aten10

R_aten10_abs <- convert_rate(rate_aten10_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1357, 
                            S = aten10$Salinity, t = aten10$Temp)
R_aten10_abs

# 28 ---------------------------------------------------------------------
inspect(aten28, time = 9, oxygen = 11, plot = TRUE)
aten28r <- data.frame(aten28$delta_t, aten28$Value)
rate_aten28 <- calc_rate(aten28r, plot=TRUE)
rate_aten28_adj <- adjust_rate(rate_aten28, bg1)

R_aten28 <- convert_rate(rate_aten28_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                        volume = 0.1328, 
                        area = 0.040998006, 
                        S = aten28$Salinity, t = aten28$Temp)
R_aten28

R_aten28_abs <- convert_rate(rate_aten28_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                            volume = 0.1328, 
                            S = aten28$Salinity, t = aten28$Temp)
R_aten28_abs

# siderastrea ---------------------------------------------------------------------
bg1 <- calc_rate.bg(sid_blank, time = 9, oxygen = 11, plot = TRUE)
bg1 # unit-less rates returned, yippee!

# 7 ---------------------------------------------------------------------
inspect(sid7, time = 9, oxygen = 11, plot = TRUE)
sid7r <- data.frame(sid7$delta_t, sid7$Value)
rate_sid7 <- calc_rate(sid7r, plot=TRUE)
rate_sid7_adj <- adjust_rate(rate_sid7, bg1)

R_sid7 <- convert_rate(rate_sid7_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                         volume = 0.1160, 
                         area = 0.020794116, 
                         S = sid7$Salinity, t = sid7$Temp)
R_sid7

R_sid7_abs <- convert_rate(rate_sid7_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                             volume = 0.1160, 
                             S = sid7$Salinity, t = sid7$Temp)
R_sid7_abs

# 15 ---------------------------------------------------------------------
inspect(sid15, time = 9, oxygen = 11, plot = TRUE)
sid15r <- data.frame(sid15$delta_t, sid15$Value)
rate_sid15 <- calc_rate(sid15r, plot=TRUE)
rate_sid15_adj <- adjust_rate(rate_sid15, bg1)

R_sid15 <- convert_rate(rate_sid15_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1137, 
                       area = 0.021581084, 
                       S = sid15$Salinity, t = sid15$Temp)
R_sid15

R_sid15_abs <- convert_rate(rate_sid15_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1137, 
                           S = sid15$Salinity, t = sid15$Temp)
R_sid15_abs

# 16 ---------------------------------------------------------------------
inspect(sid16, time = 9, oxygen = 11, plot = TRUE)
sid16r <- data.frame(sid16$delta_t, sid16$Value)
rate_sid16 <- calc_rate(sid16r, plot=TRUE)
rate_sid16_adj <- adjust_rate(rate_sid16, bg1)

R_sid16 <- convert_rate(rate_sid16_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1087, 
                       area = 0.027559014, 
                       S = sid16$Salinity, t = sid16$Temp)
R_sid16

R_sid16_abs <- convert_rate(rate_sid16_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1087, 
                           S = sid16$Salinity, t = sid16$Temp)
R_sid16_abs

# 9 ---------------------------------------------------------------------
inspect(sid9, time = 9, oxygen = 11, plot = TRUE)
sid9r <- data.frame(sid9$delta_t, sid9$Value)
rate_sid9 <- calc_rate(sid9r, plot=TRUE)
rate_sid9_adj <- adjust_rate(rate_sid9, bg1)

R_sid9 <- convert_rate(rate_sid9_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1302, 
                       area = 0.015466948, 
                       S = sid9$Salinity, t = sid9$Temp)
R_sid9

R_sid9_abs <- convert_rate(rate_sid9_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1302, 
                           S = sid9$Salinity, t = sid9$Temp)
R_sid9_abs

# 6 ---------------------------------------------------------------------
inspect(sid6, time = 9, oxygen = 11, plot = TRUE)
sid6r <- data.frame(sid6$delta_t, sid6$Value)
rate_sid6 <- calc_rate(sid6r, plot=TRUE)
rate_sid6_adj <- adjust_rate(rate_sid6, bg1)

R_sid6 <- convert_rate(rate_sid6_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1146, 
                       area = 0.017419234, 
                       S = sid6$Salinity, t = sid6$Temp)
R_sid6

R_sid6_abs <- convert_rate(rate_sid6_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1146, 
                           S = sid6$Salinity, t = sid6$Temp)
R_sid6_abs

# 12 ---------------------------------------------------------------------
inspect(sid12, time = 9, oxygen = 11, plot = TRUE)
sid12r <- data.frame(sid12$delta_t, sid12$Value)
rate_sid12 <- calc_rate(sid12r, plot=TRUE)
rate_sid12_adj <- adjust_rate(rate_sid12, bg1)

R_sid12 <- convert_rate(rate_sid12_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1231, 
                       area = 0.019507726, 
                       S = sid12$Salinity, t = sid12$Temp)
R_sid12

R_sid12_abs <- convert_rate(rate_sid12_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1231, 
                           S = sid12$Salinity, t = sid12$Temp)
R_sid12_abs

# 4 ---------------------------------------------------------------------
inspect(sid4, time = 9, oxygen = 11, plot = TRUE)
sid4r <- data.frame(sid4$delta_t, sid4$Value)
rate_sid4 <- calc_rate(sid4r, plot=TRUE)
rate_sid4_adj <- adjust_rate(rate_sid4, bg1)

R_sid4 <- convert_rate(rate_sid4_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1256, 
                       area = 0.013817342, 
                       S = sid4$Salinity, t = sid4$Temp)
R_sid4

R_sid4_abs <- convert_rate(rate_sid4_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1256, 
                           S = sid4$Salinity, t = sid4$Temp)
R_sid4_abs

# 31 ---------------------------------------------------------------------
inspect(sid31, time = 9, oxygen = 11, plot = TRUE)
sid31r <- data.frame(sid31$delta_t, sid31$Value)
rate_sid31 <- calc_rate(sid31r, plot=TRUE)
rate_sid31_adj <- adjust_rate(rate_sid31, bg1)

R_sid31 <- convert_rate(rate_sid31_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1199, 
                       area = 0.023321494, 
                       S = sid31$Salinity, t = sid31$Temp)
R_sid31

R_sid31_abs <- convert_rate(rate_sid31_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1199, 
                           S = sid31$Salinity, t = sid31$Temp)
R_sid31_abs

# 32 ---------------------------------------------------------------------
inspect(sid32, time = 9, oxygen = 11, plot = TRUE)
sid32r <- data.frame(sid32$delta_t, sid32$Value)
rate_sid32 <- calc_rate(sid32r, plot=TRUE)
rate_sid32_adj <- adjust_rate(rate_sid32, bg1)

R_sid32 <- convert_rate(rate_sid32_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1210, 
                       area = 0.022943144, 
                       S = sid32$Salinity, t = sid32$Temp)
R_sid32

R_sid32_abs <- convert_rate(rate_sid32_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1210, 
                           S = sid32$Salinity, t = sid32$Temp)
R_sid32_abs

# 25 ---------------------------------------------------------------------
inspect(sid25, time = 9, oxygen = 11, plot = TRUE)
sid25r <- data.frame(sid25$delta_t, sid25$Value)
rate_sid25 <- calc_rate(sid25r, plot=TRUE)
rate_sid25_adj <- adjust_rate(rate_sid25, bg1)

R_sid25 <- convert_rate(rate_sid25_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1217, 
                       area = 0.02368471, 
                       S = sid25$Salinity, t = sid25$Temp)
R_sid25

R_sid25_abs <- convert_rate(rate_sid25_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1217, 
                           S = sid25$Salinity, t = sid25$Temp)
R_sid25_abs

# 14 ---------------------------------------------------------------------
inspect(sid14, time = 9, oxygen = 11, plot = TRUE)
sid14r <- data.frame(sid14$delta_t, sid14$Value)
rate_sid14 <- calc_rate(sid14r, plot=TRUE)
rate_sid14_adj <- adjust_rate(rate_sid14, bg1)

R_sid14 <- convert_rate(rate_sid14_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1208, 
                       area = 0.0348082, 
                       S = sid14$Salinity, t = sid14$Temp)
R_sid14

R_sid14_abs <- convert_rate(rate_sid14_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1208, 
                           S = sid14$Salinity, t = sid14$Temp)
R_sid14_abs

# 23 ---------------------------------------------------------------------
inspect(sid23, time = 9, oxygen = 11, plot = TRUE)
sid23r <- data.frame(sid23$delta_t, sid23$Value)
rate_sid23 <- calc_rate(sid23r, plot=TRUE)
rate_sid23_adj <- adjust_rate(rate_sid23, bg1)

R_sid23 <- convert_rate(rate_sid23_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1283, 
                       area = 0.0204309, 
                       S = sid23$Salinity, t = sid23$Temp)
R_sid23

R_sid23_abs <- convert_rate(rate_sid23_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1283, 
                           S = sid23$Salinity, t = sid23$Temp)
R_sid23_abs

# 29 ---------------------------------------------------------------------
inspect(sid29, time = 9, oxygen = 11, plot = TRUE)
sid29r <- data.frame(sid29$delta_t, sid29$Value)
rate_sid29 <- calc_rate(sid29r, plot=TRUE)
rate_sid29_adj <- adjust_rate(rate_sid29, bg1)

R_sid29 <- convert_rate(rate_sid29_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1228, 
                       area = 0.03450552, 
                       S = sid29$Salinity, t = sid29$Temp)
R_sid29

R_sid29_abs <- convert_rate(rate_sid29_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1228, 
                           S = sid29$Salinity, t = sid29$Temp)
R_sid29_abs

# 27 ---------------------------------------------------------------------
inspect(sid27, time = 9, oxygen = 11, plot = TRUE)
sid27r <- data.frame(sid27$delta_t, sid27$Value)
rate_sid27 <- calc_rate(sid27r, plot=TRUE)
rate_sid27_adj <- adjust_rate(rate_sid27, bg1)

R_sid27 <- convert_rate(rate_sid27_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1200, 
                       area = 0.02845192, 
                       S = sid27$Salinity, t = sid27$Temp)
R_sid27

R_sid27_abs <- convert_rate(rate_sid27_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1200, 
                           S = sid27$Salinity, t = sid27$Temp)
R_sid27_abs

# 19 ---------------------------------------------------------------------
inspect(sid19, time = 9, oxygen = 11, plot = TRUE)
sid19r <- data.frame(sid19$delta_t, sid19$Value)
rate_sid19 <- calc_rate(sid19r, plot=TRUE)
rate_sid19_adj <- adjust_rate(rate_sid19, bg1)

R_sid19 <- convert_rate(rate_sid19_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1222, 
                       area = 0.022413454, 
                       S = sid19$Salinity, t = sid19$Temp)
R_sid19

R_sid19_abs <- convert_rate(rate_sid19_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1222, 
                           S = sid19$Salinity, t = sid19$Temp)
R_sid19_abs

# 20 ---------------------------------------------------------------------
inspect(sid20, time = 9, oxygen = 11, plot = TRUE)
sid20r <- data.frame(sid20$delta_t, sid20$Value)
rate_sid20 <- calc_rate(sid20r, plot=TRUE)
rate_sid20_adj <- adjust_rate(rate_sid20, bg1)

R_sid20 <- convert_rate(rate_sid20_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'umol/hr/cm2',
                       volume = 0.1192, 
                       area = 0.02171729, 
                       S = sid20$Salinity, t = sid20$Temp)
R_sid20

R_sid20_abs <- convert_rate(rate_sid20_adj, oxy.unit = 'hPa', time.unit = 'm',  output.unit = 'mg/min',
                           volume = 0.1192, 
                           S = sid20$Salinity, t = sid20$Temp)
R_sid20_abs


