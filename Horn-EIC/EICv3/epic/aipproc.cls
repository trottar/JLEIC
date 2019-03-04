%%%%  Start of file aipproc.cls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
\ifx\every@math@size\undefined\else              
  \let\old@expast\@expast
  \def\@expast#1{\old@expast{#1}\let\@tempa\reserved@a}
\fi

\input{aipproc.sty}
%%%%  End of file aipproc.cls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%