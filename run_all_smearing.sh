muonE=$1
echo "Smearing 0"
python make_eloss_dist.py $1 0.00
echo "Smearing 0.02"
python make_eloss_dist.py $1 0.02
echo "Smearing 0.04"
python make_eloss_dist.py $1 0.04
echo "Smearing 0.06"
python make_eloss_dist.py $1 0.06
echo "Smearing 0.08"
python make_eloss_dist.py $1 0.08
echo "Smearing 0.10"
python make_eloss_dist.py $1 0.10
echo "Smearing 0.12"
python make_eloss_dist.py $1 0.12
echo "Smearing 0.14"
python make_eloss_dist.py $1 0.14
echo "Smearing 0.16"
python make_eloss_dist.py $1 0.16
echo "Smearing 0.18"
python make_eloss_dist.py $1 0.18
echo "Smearing 0.2"
python make_eloss_dist.py $1 0.20
echo "Smearing 0.22"
python make_eloss_dist.py $1 0.22
echo "Smearing 0.24"
python make_eloss_dist.py $1 0.24
echo "Smearing 0.26"
python make_eloss_dist.py $1 0.26
echo "Smearing 0.28"
python make_eloss_dist.py $1 0.28
echo "Smearing 0.30"
python make_eloss_dist.py $1 0.30
