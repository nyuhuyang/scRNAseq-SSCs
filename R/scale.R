#2018-10-02
# scale
df <- t(data.frame("X1" = 1:5,
                   "X2" = 6:10,
                   "X3" = c(0,7,8,9, 100),
                   "X4" = c(rep(0,4),100)))
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
scale.1 <- function(df, fun1, fun2){
        df1 <- sweep(df, 1, apply(df,1, fun1),"/")
        return(sweep(df1, 1, apply(df1,1,fun2),"-"))
}
scale.1(df, fun1 = mad, fun2 = min) %>% kable() %>% kable_styling()
scale.1(df, fun1 = function(x){max(x) - min(x)}, fun2 = min) %>% kable() %>% kable_styling()
