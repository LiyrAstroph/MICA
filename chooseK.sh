for((i=8; i<=18;i++))
do
#sed -i "s/nc=.*/nc=$i;/" mcmc.c
#make
sed "s/nc.*/nc      $i/" param.txt_bak > param.txt
./rmpatch param.txt
cp data/mcmc.txt data/mcmc_$i.txt
sed -n "/aicc/p" data/results.txt | awk '{print $2 " " $3}' >> data/aicc.txt
#sed -n "/chi2/p" data/results.txt | awk '{print $2 " " $3}' >> data/chi2.txt
cd analysis/
python lcplot.py
cp fig.pdf fig_${i}.pdf
cd ../data/
cp transfer.txt transfer_${i}.txt
cd ..
done
