#2018-10-02
# scale
df <- t(data.frame("X1" = 1:5,
                   "X2" = 6:10,
                   "X3" = c(0,7,8,9, 100),
                   "X4" = c(rep(0,4),100),
                   "X5" = c(0,rep(100,4))))
library(dplyr)
library(kableExtra)
df %>% kable() %>% kable_styling()

#sights::normRobZ
df %>% t() %>% sights::normRobZ(dataRows = NULL, dataCols = NULL) %>% 
        t() %>% kable() %>% kable_styling()

df %>% t() %>% scale(center = T, scale = T) %>% t() %>% kable() %>% kable_styling()

# quantable::robustscale
df %>% quantable::robustscale(dim = 1,center = T, 
                              scale = F, 
                              preserveScale = F) -> result 
result$data %>% kable() %>% kable_styling()

# clusterSim::data.Normalization
df %>% clusterSim::data.Normalization(type="n12a",normalization="row") %>% 
        kable() %>% kable_styling()


# home made scale
.scale <- function(df, fun1 = sd, fun2 = mean){
        df1 <- sweep(df, 2, apply(df,2, fun1),"/")
        return(sweep(df1, 2, apply(df1,2,fun2),"-"))
}

df %>% t() %>% .scale(fun1 = sd, fun2 = median) %>% t() %>% kable() %>% kable_styling()

scale.1(df, fun1 = mad, fun2 = min) %>% kable() %>% kable_styling()
scale.1(df, fun1 = function(x){max(x) - min(x)}, fun2 = min) %>% kable() %>% kable_styling()
