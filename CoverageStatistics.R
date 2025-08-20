#Commands for Calculations for: Kruskal Wallis and Dunn's Test - Coverage
#Libraries ----
library(rstatix)
library(FSA)

#Coverage Statistics: ----
#Coverage Data - Median of Controls
#Protocol = P
data_BRCA2<-data.frame(
  Protocol = c(
    rep("P1", 10),
    rep("P2", 10),
    rep("P3", 10),
    rep("P4", 10),
    rep("P5", 10),
    rep("P6", 10),
    rep("P7", 2),
    rep("P8", 10)
  ),
  Coverage=c(
    41.98,65.75,70.33,47.86,68.99,77.94,57.31,70.34,53.27,76.94,
    13.62,14.11,11.54,13.05,12.99,19.59,19.24,13.19,4.71,5.78,
    16.46,28.46,26.24,14.61,40.03,20.98,89.34,39.14,44.56,54.91,
    24016.80,30919.10,17803.40,21585.40,20685.80,13299.50,19384.40,13024.00,26662.80,27946.40,
    6254.06,6240.18,6914.06,6531.28,6598.14,5960.34,6744.59,6934.77,7182.01,7358.05,
    1530.08,1165.17,1332.74,3221.22,1421.69,1174.82,1054.39,4584.47,13093.50,5583.94,
    145.27,124.23,
    37.37,44.60,57.41,67.90,31.60,57.04,112.61,84.68,59.28,84.81
  )
)

data_BRCA1<-data.frame(
  Protocol = c(
    rep("P1", 10),
    rep("P2", 10),
    rep("P3", 10),
    rep("P4", 10),
    rep("P5", 10),
    rep("P6", 10),
    rep("P7", 2),
    rep("P8", 10)
  ),
  Coverage=c(
    154.80,194.41,251.15,146.34,195.45,292.25,213.23,235.08,187.53,257.60,
    77.13,44.90,58.89,43.45,54.85,99.11,91.37,56.98,32.65,26.31,
    57.83,52.54,66.97,42.92,90.59,58.18,112.16,97.85,120.21,96.10,
    44526.10,196.04,835.30,27519.30,53445.30,44075.10,23005.40,667.49,37630.80,33398.10,
    4020.08,3308.39,3245.74,3870.37,3352.71,3861.98,3764.09,3882.39,4863.65,2940.89,
    3381.79,1752.89,3586.46,5670.82,3778.59,2662.45,1471.28,11009.80,16259.20,14112.80,
    443.11,426.59,
    53.15,79.17,57.39,101.22,34.87,59.49,114.15,88.80,62.41,83.51
  ) 
)

#Kruskal Wallis Test ----
# Perform Kruskal-Wallis test
kruskal_result_BRCA1 <- kruskal.test(Coverage ~ Protocol, data = data_BRCA1)

# Display results
print(kruskal_result_BRCA1)
#Result: Kruskal-Wallis chi-squared = 60.94, df = 7, p-value = 9.795e-11

kruskal_result_BRCA2 <- kruskal.test(Coverage ~ Protocol, data = data_BRCA2)

# Display results
print(kruskal_result_BRCA2)
#Result: Kruskal-Wallis chi-squared = 65.978, df = 7, p-value = 9.561e-12
#Conclusion: Need Dunn's Test for both BRCA1 and BRCA2

#Dunn's Test ----
dunn_BRCA1 <- dunnTest(Coverage ~ Protocol, data = data_BRCA1, method = "holm")
print(dunn_BRCA1)

#   Comparison           Z      P.unadj        P.adj
#1     P1 - P2  2.62834914 8.580040e-03 1.544407e-01
#2     P1 - P3  1.90181360 5.719553e-02 7.435419e-01
#3     P2 - P3 -0.72653553 4.675105e-01 1.000000e+00
#4     P1 - P4 -2.64971783 8.055902e-03 1.530621e-01
#5     P2 - P4 -5.27806696 1.305538e-07 3.655505e-06
#6     P3 - P4 -4.55153143 5.325684e-06 1.278164e-04
#7     P1 - P5 -2.05139445 4.022855e-02 6.436568e-01
#8     P2 - P5 -4.67974358 2.872339e-06 7.180848e-05
#9     P3 - P5 -3.95320805 7.711033e-05 1.542207e-03
#10    P4 - P5  0.59832338 5.496242e-01 1.000000e+00
#11    P1 - P6 -2.11550052 3.438731e-02 5.845843e-01
#12    P2 - P6 -4.74384966 2.096945e-06 5.661752e-05
#13    P3 - P6 -4.01731413 5.886523e-05 1.236170e-03
#14    P4 - P6  0.53421730 5.931912e-01 1.000000e+00
#15    P5 - P6 -0.06410608 9.488858e-01 9.488858e-01
#16    P1 - P7 -0.40095966 6.884498e-01 1.000000e+00
#17    P2 - P7 -1.91843774 5.505553e-02 7.707774e-01
#18    P3 - P7 -1.49897225 1.338808e-01 1.000000e+00
#19    P4 - P7  1.12885565 2.589587e-01 1.000000e+00
#20    P5 - P7  0.78341348 4.333843e-01 1.000000e+00
#21    P6 - P7  0.82042514 4.119738e-01 1.000000e+00
#22    P1 - P8  2.04071010 4.127965e-02 6.191948e-01
#23    P2 - P8 -0.58763903 5.567746e-01 1.000000e+00
#24    P3 - P8  0.13889650 8.895319e-01 1.000000e+00
#25    P4 - P8  4.69042793 2.726342e-06 7.088490e-05
#26    P5 - P8  4.09210455 4.274758e-05 9.404468e-04
#27    P6 - P8  4.15621062 3.235695e-05 7.442100e-04
#28    P7 - P8  1.57916418 1.142984e-01 1.000000e+00

dunn_BRCA2 <- dunnTest(Coverage ~ Protocol, data = data_BRCA2, method = "holm")
print(dunn_BRCA2)

#Comparison          Z      P.unadj        P.adj
#1     P1 - P2  2.5108213 1.204506e-02 2.168112e-01
#2     P1 - P3  1.1859624 2.356371e-01 1.000000e+00
#3     P2 - P3 -1.3248589 1.852180e-01 1.000000e+00
#4     P1 - P4 -4.0600515 4.906189e-05 1.079362e-03
#5     P2 - P4 -6.5708728 5.002117e-11 1.400593e-09
#6     P3 - P4 -5.2460139 1.554252e-07 4.041054e-06
#7     P1 - P5 -2.8954578 3.786059e-03 7.193513e-02
#8     P2 - P5 -5.4062791 6.434750e-08 1.737382e-06
#9     P3 - P5 -4.0814202 4.476135e-05 1.029511e-03
#10    P4 - P5  1.1645937 2.441835e-01 1.000000e+00
#11    P1 - P6 -2.0513944 4.022855e-02 5.631997e-01
#12    P2 - P6 -4.5622158 5.061659e-06 1.265415e-04
#13    P3 - P6 -3.2373569 1.206424e-03 2.533491e-02
#14    P4 - P6  2.0086571 4.457352e-02 5.794557e-01
#15    P5 - P6  0.8440633 3.986340e-01 1.000000e+00
#16    P1 - P7 -0.7464018 4.554247e-01 1.000000e+00
#17    P2 - P7 -2.1960252 2.809014e-02 4.494422e-01
#18    P3 - P7 -1.4311175 1.523965e-01 1.000000e+00
#19    P4 - P7  1.5976700 1.101164e-01 1.000000e+00
#20    P5 - P7  0.9252915 3.548143e-01 1.000000e+00
#21    P6 - P7  0.4379713 6.614071e-01 1.000000e+00
#22    P1 - P8  0.1068435 9.149132e-01 9.149132e-01
#23    P2 - P8 -2.4039779 1.621776e-02 2.757018e-01
#24    P3 - P8 -1.0791190 2.805347e-01 1.000000e+00
#25    P4 - P8  4.1668950 3.087767e-05 7.410640e-04
#26    P5 - P8  3.0023012 2.679469e-03 5.358938e-02
#27    P6 - P8  2.1582379 3.090934e-02 4.636401e-01
#28    P7 - P8  0.8080879 4.190400e-01 1.000000e+00

# BRCA1 significant differences
cat("=== BRCA1 SIGNIFICANT DIFFERENCES (P.adj < 0.05) ===\n")
brca1_significant <- dunn_BRCA1$res[dunn_BRCA1$res$P.adj < 0.05, ]
if(nrow(brca1_significant) > 0) {
  print(brca1_significant)
  cat("\nNumber of significant comparisons for BRCA1:", nrow(brca1_significant), "\n")
} else {
  cat("No significant differences found for BRCA1\n")
}

cat("\n")

# BRCA2 significant differences
cat("=== BRCA2 SIGNIFICANT DIFFERENCES (P.adj < 0.05) ===\n")
brca2_significant <- dunn_BRCA2$res[dunn_BRCA2$res$P.adj < 0.05, ]
if(nrow(brca2_significant) > 0) {
  print(brca2_significant)
  cat("\nNumber of significant comparisons for BRCA2:", nrow(brca2_significant), "\n")
} else {
  cat("No significant differences found for BRCA2\n")
}

# Combined
all_significant <- rbind(
  data.frame(Gene = "BRCA1", brca1_significant),
  data.frame(Gene = "BRCA2", brca2_significant)
)

cat("\n=== COMBINED SIGNIFICANT RESULTS ===\n")
print(all_significant)

#Result: Significant differences: 
#=== COMBINED SIGNIFICANT RESULTS ===
# Gene Comparison         Z      P.unadj        P.adj
#5   BRCA1    P2 - P4 -5.278067 1.305538e-07 3.655505e-06
#6   BRCA1    P3 - P4 -4.551531 5.325684e-06 1.278164e-04
#8   BRCA1    P2 - P5 -4.679744 2.872339e-06 7.180848e-05
#9   BRCA1    P3 - P5 -3.953208 7.711033e-05 1.542207e-03
#12  BRCA1    P2 - P6 -4.743850 2.096945e-06 5.661752e-05
#13  BRCA1    P3 - P6 -4.017314 5.886523e-05 1.236170e-03
#25  BRCA1    P4 - P8  4.690428 2.726342e-06 7.088490e-05
#26  BRCA1    P5 - P8  4.092105 4.274758e-05 9.404468e-04
#27  BRCA1    P6 - P8  4.156211 3.235695e-05 7.442100e-04
#4   BRCA2    P1 - P4 -4.060052 4.906189e-05 1.079362e-03
#51  BRCA2    P2 - P4 -6.570873 5.002117e-11 1.400593e-09
#61  BRCA2    P3 - P4 -5.246014 1.554252e-07 4.041054e-06
#81  BRCA2    P2 - P5 -5.406279 6.434750e-08 1.737382e-06
#91  BRCA2    P3 - P5 -4.081420 4.476135e-05 1.029511e-03
#121 BRCA2    P2 - P6 -4.562216 5.061659e-06 1.265415e-04
#131 BRCA2    P3 - P6 -3.237357 1.206424e-03 2.533491e-02
#251 BRCA2    P4 - P8  4.166895 3.087767e-05 7.410640e-04


