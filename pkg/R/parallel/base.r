#!/usr/bin/env Rscript

library(foreach)

data(geneData, package = 'Biobase')

pair <- combn(1:nrow(geneData), 2, simplify = F)

fakeData <- cbind(geneData, geneData, geneData, geneData)

pair2 <- sample(pair, 300)

print(system.time(

                  out <- foreach(p = pair2, .packages = 'boot', .combine = 'rbind') %dopar%

                  {

                    mydata <- cbind(fakeData[p[1],], fakeData[p[2], ])

                    mycor <- function(x, i) cor(x[i,1], x[i,2])

                    boot.out <- boot(mydata, mycor, 1000)

                    ci <- boot.ci(boot.out, type = 'bca')$bca[4:5]

                    c(p, ci)

                  }

                  ))

print(head(out)) # print the head of the result
