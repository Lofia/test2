**Run the following code in R console to install:**
if(!require("devtools")){install.packages("devtools");library(devtools)}
install_github('Lofia/test2')

**Try the following code:**
library(test2)
?power_of_P_test
power_of_P_test('ux',c(10,20,30),c(0.2,0.3))
?test1
test1(seq(0,1,0.01),alpha=0.05)
