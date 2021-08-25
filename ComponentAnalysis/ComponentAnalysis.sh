R="GM12878_1"
Q="GM12878_2"
O="GM12878_1vs2"
PYTHONVERSION="python3.8"

# backup
cp ~/.local/lib/$PYTHONVERSION/site-packages/skimage/metrics/_structural_similarity.py ~/backup.py

# luminance
cp ssim_luminance.py ~/.local/lib/$PYTHONVERSION/site-packages/skimage/metrics/_structural_similarity.py
chess sim -p 4 ../Processed/$R.chr2.25kb.KR.cool ../Processed/$Q.chr2.25kb.KR.cool ../Pairs/hg19.chr2.2Mb.500kb.txt ../CHESS/chr2.2Mb.500kb/$O.luminance.tsv

# contrast
cp ssim_contrast.py ~/.local/lib/$PYTHONVERSION/site-packages/skimage/metrics/_structural_similarity.py
chess sim -p 4 ../Processed/$R.chr2.25kb.KR.cool ../Processed/$Q.chr2.25kb.KR.cool ../Pairs/hg19.chr2.2Mb.500kb.txt ../CHESS/chr2.2Mb.500kb/$O.contrast.tsv

# structure
cp ssim_structure.py ~/.local/lib/$PYTHONVERSION/site-packages/skimage/metrics/_structural_similarity.py
chess sim -p 4 ../Processed/$R.chr2.25kb.KR.cool ../Processed/$Q.chr2.25kb.KR.cool ../Pairs/hg19.chr2.2Mb.500kb.txt ../CHESS/chr2.2Mb.500kb/$O.structure.tsv

# recovery
cp ~/backup.py ~/.local/lib/$PYTHONVERSION/site-packages/skimage/metrics/_structural_similarity.py

