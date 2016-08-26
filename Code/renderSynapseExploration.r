require(rmarkdown)
system.time(render("SynapseExploration.Rmd"))
system('say -r 200 Your R-script has completed')
