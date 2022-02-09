R="NORMAL"
Q="DLBCL"
O="DLBCLvsNORMAL"
PYTHONVERSION="python3.7"

# backup
cp ~/miniconda3/lib/$PYTHONVERSION/site-packages/chess/sim.py sim_backup.py

# inverse of Fano factor ("SN")
chess sim -p 4 ../COOL/$R.chr2.25kb.KR.cool ../COOL/$Q.chr2.25kb.KR.cool ../Pairs/hg19.chr2.2Mb.500kb.txt ../CHESS/chr2.2Mb.500kb/$O.tsv

# foldchange
cp sim_foldchange.py ~/miniconda3/lib/$PYTHONVERSION/site-packages/chess/sim.py
chess sim -p 4 ../COOL/$R.chr2.25kb.KR.cool ../COOL/$Q.chr2.25kb.KR.cool ../Pairs/hg19.chr2.2Mb.500kb.txt ../CHESS/chr2.2Mb.500kb/$O.foldchange.tsv

# variance
cp sim_variance.py ~/miniconda3/lib/$PYTHONVERSION/site-packages/chess/sim.py
chess sim -p 4 ../COOL/$R.chr2.25kb.KR.cool ../COOL/$Q.chr2.25kb.KR.cool ../Pairs/hg19.chr2.2Mb.500kb.txt ../CHESS/chr2.2Mb.500kb/$O.variance.tsv

# recovery
cp sim_backup.py ~/miniconda3/lib/$PYTHONVERSION/site-packages/chess/sim.py

