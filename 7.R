library('Biostrings')

Score = function(s, DNA, l){
  block = DNAString()
  for (seq_idx in 1:length(DNA)){
    block = c(block, subseq(DNA[seq_idx], start=s[seq_idx], width=l))
  }
  frequency = consensusMatrix(block)
  score = 0
  for (seq_idx in 1:ncol(frequency)){
    score = score + max(frequency[,seq_idx])
  }
}


NextLeaf = function(a, L, k){
  for (i in L:1){
    if (a[i] < k){
      a[i] = a[i] + 1
      return(a)
    }
    a[i] = 1
  }
  return(a)
}


BFMotifSearch = function(DNA, t, n, l){
  bestMotif = c()
  s = c(1)
  bestScore = Score(s, DNA, l)
  while (TRUE){
    s = NextLeaf(s, t, n-l+1)
    if (Score(s, DNA, l) > bestScore){
      bestMotif = append(bestMotif, s)
    }
    if (unique(s) == 1) return(bestMotif)
  }
}


NextVertex = function(a, i , L, k){
  if (i < L){
    a[i+1] = 1
    return(c(a, i+1))
  } 
  else{
    for (j in L:1){
      if (a[j] < k){
        a[j] = a[j] + 1
        return(c(a,j))
      }
    }
  }
  return(c(a,0))
}


ByPass = function(a, i, L, k){
  for (j in i:1){
    if (a[j] < k){
      a[j] = a[j] + 1
      return(c(a,j))
    }
  }
  return(c(a,0))
}


BBMotifSearch = function(DNA, t, n, l){
  s = (1)
  bestmotif = c()
  bestScore = 0
  i = 1
  while (i > 0){
    if (i < t){
      optimisticScore = Score(s, i, DNA, l) + (t-i)*l
      if (optimisticScore < bestScore){
        x = ByPass(s, i, t, n-l+1)
        s = x[1]
        i = x[2]
      } 
      else{
        x = NextVertex(s, i, t, n-l+1)
        s = x[1]
        i = x[2]
      }
    }
    else {
      if(Score(s, t, DNA, l) > bestScore){
      bestScore = Score(s, t, DNA, l)
      bestMotif = append(bestMotif, s[i])
      }
      x = NextVertex(s, i, t, n-l+1)
      s = x[1]
      i = x[2]
    }
  }
  return(bestmotif)
}
