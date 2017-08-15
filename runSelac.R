rm(list=ls())


## Make sure SELAC is installed with a version of R that has the 
## appropriate BLAS libraries.
## The version on Gauley that I've been using is found in the directory /opt/R/3.3.1-2.
library(selac)


## Function for reading an input file 
## File must have a structure like the example below.

## Example of input file structure. 

## Run,Working_Directory,Tree_File,Codon_Data_Path,Output_Path,Num_Cores_Gene,Num_Cores_Site,Optimizer
## 1,./,tree_rokas.tre,aln_rokas_cut/,aln_rokas_cut/SBPLX/,4,1,NLOPT_LN_SBPLX
## 2,./,tree_rokas.tre,aln_rokas_cut/,aln_rokas_cut/PRAXIS/,4,1,NLOPT_LN_PRAXIS

## All values are read in as characters, as opposed to allowing R to convert them to what
## it determines to be the most appropriate type. This prevents the list from being
## converted to factors when passed to runSelac function

## If you decide to add more columns in the input file, this function will need to be modified.
readInput <- function(input.file)
{
  runs <- read.table(input.file,sep=",",header=T,colClasses=c("character","character","character","character","character","character","character","character"))
  return(runs)
}

## Function for visually comparing results for various optimizers
## Parameters:
## 1) results: row bound list of results from SELAC runs using different optimizers 
## 2) what: vector of comparisons to makes
## 3) pdf.out: create pdf
## 4) pdf.name: name of pdf file to be output
## 5) ...: plotting parameters, such as cex.axis
plotOptimizerComparisons <- function(results,what=c("Runtime","LogLikelihood","AIC"),pdf.out = TRUE, pdf.name ="optimizer_comparisons.pdf",...)
{
  if(pdf.out){
    pdf(pdf.name)
  }
  opt.vector<-unlist(lapply(results[,"opts"],function(j){j$algorithm}))
  #cat(factor(opt.vector),"\n")
  if("Runtime" %in% what)
  {
    plot(factor(opt.vector),unlist(results[,"runtimes"]),xlab="",ylab="Runtimes (hours)",main="Runtime",...)
  }
  if("LogLikelihood" %in% what)
  {
    plot(factor(opt.vector),unlist(results[,"loglik"]),xlab="",ylab="LogLikelihood",main="LogLikelihood",...)
  }
  if("AIC" %in% what)
  {
    plot(factor(opt.vector),unlist(results[,"AIC"]),xlab="",ylab="AIC",main="AIC",...)
  }
  if(pdf.out)
  {
    dev.off()
  }

}

## Function for visually comparing results for different number of cores at either per
## gene or per site
## Parameters:
## 1) results: row bound list of results from SELAC runs using different numbers of cores
## 2) per.gene: if TRUE, compare runtimes for different number of cores per gene; else, compare runtimes
##              for different number of cores per site
## 3) pdf.out: create pdf
## 4) pdf.name: name of pdf file to be output
## 5) ...: plotting parameters, such as main and cex.axis
plotRuntimeComparisonForDiffCores <- function(results,per.gene=TRUE,pdf.out=TRUE,pdf.name="runtime_comparisons.pdf",...)
{
  if(per.gene){
    ncores <- unlist(results[,"ncores.per.gene"])
    xlab <- "Number of Cores Per Gene"
  } else{
    ncores <- unlist(results[,"ncores.per.site"])
    xlab <- "Number of Cores Per Site"
  }
  plot(ncores,unlist(results.data[,"runtimes"]),xlab=xlab,ylab="Runtimes (hours)")

}

## Function for creating a row-bound list of saved SELAC outputs
## Parameters:
## 1) rds.files: vector of paths to saved rds files from runSELAC
extractOutputFromRDSFiles <- function(rds.files)
{
  output <- lapply(rds.files,function(file.i){results<-readRDS(file.i)})
  output<- do.call("rbind",output)
  return(output)
}


## Runs SELAC
## Parameters:
## 1) arg.list: list of arguments either provided by the user via command line or 
##              with an input file
runSelac <- function(arg.list)
{
  run.number <- as.numeric(arg.list[1])
  directory <- arg.list[2]
  tree.file <- arg.list[3]
  codon.data.path <- arg.list[4]
  out.path <- arg.list[5]
  ncores.per.gene <- as.numeric(arg.list[6])
  ncores.per.site <- as.numeric(arg.list[7])
  opt.method <-arg.list[8]
  
  cat("\nRun Number: ",run.number,"\nOptimization Algorithm: ",opt.method)
  setwd(directory)
  set.seed(123456)
  # read tree
  tree <- read.tree(tree.file)
  
  cat("\n\n ### ncores.per.gene: ", ncores.per.gene, "\t ncores.per.site: ", ncores.per.site, "\n\n")
  # run selac
  cat("Setup Done. Starting SELAC\n")
  start <- Sys.time()
  result <- SelacOptimize(codon.data.path, n.partitions=5, phy=tree, data.type="codon", codon.model="selac", edge.length="optimize", 
                          edge.linked=TRUE, optimal.aa="optimize", nuc.model="GTR", include.gamma=TRUE, gamma.type="quadrature", 
                          ncats=4, numcode=1, diploid=TRUE, k.levels=0, aa.properties=NULL, verbose=FALSE, n.cores.by.gene=ncores.per.gene, 
                          n.cores.by.gene.by.site=ncores.per.site, max.tol=1e-4, max.evals=10, max.restarts=1, 
                          user.optimal.aa=NULL, fasta.rows.to.keep=NULL, recalculate.starting.brlen=TRUE, 
                          output.by.restart=TRUE, output.restart.filename="aln_rokas_cut", 
                          user.supplied.starting.param.vals=NULL, tol.step=1,optimizer.algorithm=opt.method)
  end <- Sys.time()
  runtime <- end - start
  cat("\n\n ### Runtime: ", runtime, "\n\n")
  result <- c(result,runtimes=as.numeric(runtime,units="hours"),ncores.per.gene=ncores.per.gene,ncores.per.site=ncores.per.site)
  
  # save the results of SELAC to a file so can easily be later accessed via the
  # extractOutputFromRDSFiles function
  saveRDS(result,paste0(out.path,"selac_out.rds"))  
                                                  
  return(result)
}

## Can call this script via the command:
## Rscript runSelac [arguments]
## Examples:

## 1) With just input file:
##    Rscript runSelac.R "input.file='input.txt'"

## 2) With input file and run in parallel:
##    Rscript runSelac.R "input.file='input.txt'" "run.parallel=TRUE"

## 3) 
##    Rscript runSelac.R "directory='./'" "tree.file='tree_rokas.tre'" "codon.data.path='aln_rokas_cut/'" "out.path='aln_rokas_cut/'" "ncores.per.gene=4" 
##           "ncores.per.site=1" "opt.method='NL_LN_SBPLX'"

args<-(commandArgs(TRUE))
input.file <- NULL
if(length(args)==0) {
  # Note: if you want to run this script via Rstudio, you can comment out the stop call and initialize
  # the appropriate variables
  stop("At least an input file must be provided to runSelac\n",call.=T)
  #input.file <- "input.txt"
  #run.parallel <- FALSE
} else if(length(args)==1) {
  eval(parse(text=args[[1]]))
  if(is.null(input.file))
  {
    stop("Must provide at minimum an input file to run SELAC\n",call.=T)
  }
  run.parallel <- FALSE
} else if(length(args)==2){
  for(i in 1:length(args))
  {
    eval(parse(text=args[[i]]))
  }
} else if(length(args)==7){
  for(i in 1:length(args))
  {
    eval(parse(text=args[[i]]))
  }
  run.number <- "1"
} else {
  stop("You are missing 1 or more command line arguments. Check that you have provided values for:\n1)directory
  \n2)tree.file\n3)codon.data.path\n4)out.path\n5)ncores.per.gene\n6)ncores.per.site\n7)opt.method",call.=FALSE)
}


if (!is.null(input.file))
{
  runs <- readInput(input.file)

  ## Transpose the data frame because lapply works on columns, not rows. Input file uses rows because I felt
  ## they would be easier to read and modify.
  char <- rep("character",length(nrow(runs)))
  runs.t <- data.frame(t(runs),stringsAsFactors = FALSE)
  if(run.parallel)
  {
    library(parallel)
    ## User will have to be aware of the number of cores on machine. My understanding is 
    ## mclappy and similar functions cannot communicate between nodes on a machine, only 
    ## between cores on the same node. Potentially, we will need to look into libraries
    ## allowing R to interface with MPI libraries. 

    ## If I'm understanding how SELAC and mclapply work, a process will be created when this file 
    ## is sourced and then n child processes will be created, one for each row in your input
    ## file. Each of these child processes will also fork their own child processes based on the number
    ## of cores per gene and per site. This could be come problematic if trying to run a large number of 
    ## jobs in parallel on a computer or node with relatively few cores. For example, Beacon compute nodes each have 16
    ## cores, so a run with 4 jobs (4 rows in the input files) with 5 cores per gene and 1 core per site will
    ## take up more cores than are avaiable. This will result in jobs sharing cores, causing increased overhead.

    ## Note that I have set mc.preschedule to FALSE. Taken from the R documentation on mclapply

    ## "if set to TRUE then the computation is first divided to (at most) as many jobs are there are cores and then the jobs are started, 
    ## each job possibly covering more than one value. If set to FALSE then one job is forked for each value of X. The former is better 
    ## for short computations or large number of values in X, 
    ## the latter is better for jobs that have high variance of completion time and not too many values of X compared to mc.cores."

    ## It seems to me the latter case would be more common if using a machine like Gauley. 
    results <- mclapply(runs.t,runSelac,mc.preschedule=FALSE,mc.cores=nrow(runs))
  } else{
    results <- lapply(runs.t,runSelac)
  }
  results<-do.call("rbind",results)
  plotOptimizerComparisons(results,pdf.name="test.pdf")
} else{
  arguments <- c(run.number,directory,tree.file,codon.data.path,out.path,ncores.per.gene,ncores.per.site,opt.method)
  results <- runSelac(arguments)
}
