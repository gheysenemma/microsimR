#' @title Generate an interaction matrix
#'
#' @description Generate an interaction matrix, either randomly from a uniform distribution or
#' using Klemm-Eguiluz algorithm to generate a modular and scale-free interaction matrix.
#'
#' @param N number of species
#' @param type random (sample a uniform distribution), klemm (generate a Klemm-Eguiluz matrix) or empty (zero everywhere, except for diagonal which is set to d)
#' @param pep desired positive edge percentage (only for klemm)
#' @param d diagonal values (should be negative)
#' @param min.strength random: minimal off-diagonal interaction strength (only for random)
#' @param max.strength random: maximal off-diagonal interaction strength, klemm: maximal absolute off-diagonal interaction strength
#' @param c desired connectance (interaction probability)
#' @param ignore.c do not adjust connectance
#' @param negedge.symm set symmetric negative interactions (only for klemm)
#' @param clique.size modularity parameter (only for klemm)
#' @param groups vector of group memberships for each species, assign NA if species does not belong to any group (only for random)
#' @param intra.group.strength interaction strength between members of the same group
#' @param inter.group.strength interaction strength between members of different groups (if not defined, will be assigned randomly)
#' @return the interaction matrix
#' @examples
#' Adefault = randomA(N=100, d=-0.5, min.strength=-0.5, max.strength=0.5, c=0.02, groups=c(), intra.group.strength=0.5, inter.group.strength=NA)
#' groups=c(rep(NA,5),rep(1,10),rep(2,5),rep(3,10),rep(4,10))
#' Agroup=randomA(N=40,groups=groups,c=0.5,intra.group.strength=0.1,inter.group.strength=-0.5, d=-1)
#' @references Klemm & Eguiluz, Growing Scale-Free Networks with Small World Behavior \url{http://arxiv.org/pdf/cond-mat/0107607v1.pdf}
#' @note TEMPORARILY MERGE OF COPIES FROM SEQTIME FUNCTION, NOT YET OPTIMISED
#' @export

randomA<-function(N=100, type="random",pep=50, d=-0.5, min.strength=-0.5, max.strength=0.5, c=0.02, ignore.c=FALSE, negedge.symm=FALSE, clique.size=5, groups=c(), intra.group.strength=0.5, inter.group.strength=NA){
  A=matrix(0,nrow=N,ncol=N)      # init species interaction matrix
  if(type=="random"){
    if(length(groups)>0){
      if(length(groups)!=nrow(A)){
        stop("Please define a group membership for each species.")
      }
    }
    for (i in 1:N){
      for(j in 1:N){
        if(i==j){
          A[i,j]=d
        }else{
          if(length(groups)==0){
            A[i,j] = runif(1,min=min.strength,max=max.strength)
          }else{
            group1=groups[i]
            group2=groups[j]
            if(!is.na(group1) && !is.na(group2) && group1==group2){
              A[i,j] = intra.group.strength
            }else{
              # assign interaction strength between groups randomly
              if(is.na(inter.group.strength)){
                A[i,j] = runif(1,min=min.strength,max=max.strength)
              }else{
                # assign selected interaction strength between groups
                A[i,j] = inter.group.strength
              }
            }
          }
        }
      }
    }
  }else if(type=="empty"){
    diag(A)=d
  }else if(type=="klemm"){
    g<-klemm.game(N,verb=FALSE, clique.size)
    A=get.adjacency(g)
    A=as.matrix(A)
    diag(A)=d
  }

  if(ignore.c==FALSE){
    #    print(paste("Adjusting connectance to",c))
    A=modifyA(A,c=c, mode="adjustc")
  }

  # for Klemm-Eguiluz: inroduce negative edges and set random interaction strengths
  if(type=="klemm"){
    if(pep < 100){
      A=modifyA(A=A, perc=(100-pep), symmetric=negedge.symm, mode="negpercent")
    }
    # excluding edges on the diagonal
    #    print(paste("Final arc number (excluding self-arcs)", length(A[A!=0])-N ))
    # excluding negative edges on the diagonal
    #    print(paste("Final negative arc number (excluding self-arcs)", length(A[A<0])-N ))

    # check PEP and number of asymmetric negative interactions
    # assuming diagonal values are negative
    pep = getPep(A)
    #    print(paste("PEP:",pep))

    # convert binary interaction strengths (-1/1) into continuous ones using uniform distribution
    # zero would remove the edge, so the minimum strength is small, but non-zero
    min.klemm.strength=0.00001
    for(i in 1:nrow(A)){
      for(j in 1:nrow(A)){
        # skip diagonal
        if(i != j){
          A[i,j]=A[i,j]*runif(1,min=min.klemm.strength,max=max.strength)
        }
      }
    }
  }
  return(A)
}


modifyA<-function(A, mode="adjustc", strength="binary", factor=2, minstrength=0.1, c=0.2, perc=50, symmetric=FALSE){
  edgeNumAdded = 0
  #  print(paste("Initial edge number", length(A[A!=0])))
  c_obs = getConnectance(A)
  #  print(paste("Initial connectance", c_obs))
  if(mode == "adjustc"){
    if(c_obs < c){
      while(c_obs < c){
        # randomly select source node of edge
        xpos=sample(c(1:ncol(A)))[1]
        # randomly select target node of edge
        ypos=sample(c(1:ncol(A)))[1]
        # avoid diagonal
        if(xpos != ypos){
          # count as added if there was no edge yet
          if(A[xpos,ypos]==0){
            edgeNumAdded = edgeNumAdded+1
          }
          # add edge
          A[xpos,ypos]=getStrength(strength=strength,pos=TRUE, minstrength=minstrength)
          c_obs=getConnectance(A=A)
        }
      }
      #      print(paste("Number of edges added", edgeNumAdded))
    }else if(c_obs > c){
      edgeNumRemoved = 0
      while(c_obs > c){
        xpos=sample(c(1:ncol(A)))[1]
        ypos=sample(c(1:ncol(A)))[1]
        # avoid diagonal
        if(xpos != ypos){
          # count as removed if there was an edge before
          if(A[xpos,ypos]!=0){
            edgeNumRemoved = edgeNumRemoved+1
          }
          # remove edge
          A[xpos,ypos]=0
          c_obs = getConnectance(A)
        }
      }
      #      print(paste("Number of edges removed", edgeNumRemoved))
    }
    #    print(paste("Final connectance", c_obs))
  }else if(mode=="mergeposlinks" || mode=="mergeneglinks" || mode=="mergelinks"){
    taxonnames=as.character(A[,1])
    entries=unique(taxonnames)
    # remove taxon name column
    A=A[,2:ncol(A)]
    mergedlinks=matrix(0,nrow=length(entries),ncol=length(entries))
    rownames(mergedlinks)=entries
    colnames(mergedlinks)=entries
    for(i in 1 : nrow(A)){
      for(j in 1 : ncol(A)){
        if(!is.na(A[i,j]) && A[i,j]!=0){
          merge=FALSE
          xIndex=which(entries==taxonnames[i])
          yIndex=which(entries==taxonnames[j])
          #print(paste("x index:",xIndex))
          #print(paste("y index:",yIndex))
          if(mode=="mergeposlinks" && A[i,j]>0){
            merge=TRUE
          }else if(mode=="mergeneglinks" && A[i,j]<0){
            merge=TRUE
          }else if(mode=="mergelinks"){
            merge=TRUE
          }
          if(merge==TRUE){
            if(mode=="mergelinks"){
              if(A[i,j]<0){
                mergedlinks[xIndex,yIndex]=mergedlinks[xIndex,yIndex]-1
              }else{
                mergedlinks[xIndex,yIndex]=mergedlinks[xIndex,yIndex]+1
              }
            }else{
              mergedlinks[xIndex,yIndex]=mergedlinks[xIndex,yIndex]+1
            }
          }
        } # interaction is not zero
      } # column loop
    } # row loop
    A=mergedlinks
  }else if(mode=="removeorphans"){
    # since A can be asymmetric, only those species can be removed for which rows and columns are simultaneously zero (except for diagonal)
    toKeep=c()
    diagvals=diag(A)
    diag(A)=0
    for(i in 1:nrow(A)){
      rowsum=sum(abs(A[i,]))
      colsum=sum(abs(A[,i]))
      if(rowsum != 0 || colsum!=0){
        toKeep=append(toKeep,i)
      }
    }
    A=A[toKeep,toKeep]
    diag(A)=diagvals[toKeep]
  }else if(mode == "negpercent"){
    # arc number: number of non-zero entries in the interaction matrix
    num.edge=length(A[A!=0])
    num.perc=(num.edge/100)*perc
    # subtract the negative edges that are already there
    num.neg.edge=length(A[A<0])
    #    print(paste("Number of negative edges already present:",num.neg.edge))
    if(num.neg.edge>num.perc){
      warning("The matrix has more negative edges than are required to reach the desired negative edge percentage!")
    }else if(num.neg.edge==num.perc){
      #      print("The matrix has already the desired negative edge percentage.")
    }else{
      # those negative edges already present do not need to be added
      num.perc=num.perc-num.neg.edge
      # symmetric interactions: we will count negative edges, not arcs
      if(symmetric == TRUE){
        num.perc=round(num.perc/2)
      }else{
        num.perc=round(num.perc)
      }
      #      print(paste("Converting",num.perc,"edges into negative edges",sep=" "))
      indices=which(A>0,arr.ind=T)
      # randomly select indices
      rand=sample(c(1:nrow(indices)))
      xyposAlreadySeen = c()
      counter = 0
      # loop over number of negative edges to introduce
      for(i in 1:num.perc){
        xpos=indices[rand[i],1]
        ypos=indices[rand[i],2]
        xypos=paste(xpos,"_",ypos, sep="")
        yxpos=paste(ypos,"_",xpos,sep="")
        # if we find an index pair that was already used, we have to look for another index pair,
        # since using the same index pair means to use the same arc or the same arc in reverse direction
        if(symmetric == TRUE && is.element(xypos,xyposAlreadySeen) == TRUE){
          xpos = indices[rand[nrow(indices)-counter],1]
          ypos = indices[rand[nrow(indices)-counter],2]
          counter = counter + 1
          if((num.perc + counter) > nrow(indices)){
            stop("More negative edges requested than can be set!")
          }
        }
        xyposAlreadySeen = c(xypos, yxpos, xyposAlreadySeen)
        # print for tests
        # print(paste("xpos",xpos,"ypos",ypos,"value:",A[xpos,ypos],sep=" "))
        negval=getStrength(strength=strength,pos=FALSE,minstrength=minstrength)
        A[xpos,ypos]=negval
        if(symmetric == TRUE){
          A[ypos,xpos]=negval
        }
        #print(paste("xpos",xpos,"ypos",ypos,"value:",A[xpos,ypos],sep=" "))
        #print(paste("reverse value:",A[ypos,xpos],sep=" "))
      }
    }
  }else if(mode == "tweak"){
    # check that positive arcs are present
    if(length(A[A>0]) > 0){
      # row and column indices of positive arcs
      indices.pos = which(A>0,arr.ind=TRUE)
      # randomly select a positive arc
      x=sample(1:nrow(indices.pos),1)
      # convert positive arc into negative one, keeping the same interaction strength
      A[indices.pos[x,1],indices.pos[x,2]]=A[indices.pos[x,1],indices.pos[x,2]]*(-1)
    }else{
      warning("Cannot tweak. No positive arc in the given matrix.")
    }
  }else if(mode == "enforceneg"){
    diag=diag(A)
    indices.neg = which(A<0,arr.ind=TRUE)
    # multiply negative entries by given factor
    A[indices.neg]=A[indices.neg]*factor
    # keep original diagonal
    diag(A)=diag
  }else if(mode == "schur"){
    # remove positive real parts of eigenvalues if any (using schur decomposition)
    sd<-dim(A)

    if(max(Re(eigen(A)$values))){
      # division by max.A helps removing positive eigenvalues
      max=max(A)
      A=A/max

      diagt<-diag(sd[2])+0i

      # Computes the generalized eigenvalues and Schur form of a pair of matrices.
      # R=imaginary part identical to 0 with a tolerance of 100*machine_precision as determined by Lapack
      schur.A<-geigen::gqz(A,diagt,"R")
      # generalized inverse of a matrix
      T<-schur.A$S%*%MASS::ginv(schur.A$T)
      rediag<-Re(diag(T))
      imdiag<-Im(diag(T))

      indicesP=rediag>0
      listind=1:sd[2]

      for(k in listind[indicesP]){
        T[k,k]<- complex(real=-Re(T[k,k]),imaginary=Im(T[k,k]))
      }

      A <- schur.A$Q %*% T %*% MASS::ginv(schur.A$Q)
      A<-Re(A)
      A=A*max
    }

  }else{
    stop(paste("Mode",mode,"not known."))
  }
  c=getConnectance(A)
  #  print(paste("Final connectance:",c))
  return(A)
}

getConnectance <- function(A){

  N <- nrow(A)

  # exclude diagonal from observed and possible interactions
  c <- (length(A[A!=0])-N)/(ncol(A)*ncol(A)-N)

  return(c)
}

################## helper functions ################

# Get the interaction strength.
getStrength<-function(strength="binary", minstrength=0.1, pos=TRUE){
  value = NA
  if(strength=="binary"){
    value = 1
  }else if(strength == "uniform"){
    value = runif(1,min=minstrength,max=1)
  }
  if(!pos){
    value = -1*value
  }
  return(value)
}

# Check whether the interaction matrix is fully connected (entirely filled with 1)
isFullyconnected<-function(A){
  if(length(A[A!=0])==nrow(A)*nrow(A)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

