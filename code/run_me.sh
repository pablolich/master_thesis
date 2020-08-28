#First, simulate communities with cross feeding interactions
echo -e "\n1. COMMUNITY ASSEMBLY"
echo -e "\n Running assembly simulations"
python3 community_assembly.py 
python3 interaction_evolution.py
echo -e "\n Analysis of simulated communities"
Rscript analysis_community_assembly.r 
echo -e "\n2. COMMUNITY COALESCENCE"
echo -e "\n Running coalescencee simulations"
python3 tile_plot.py
python3 s_plot.py 
echo -e "\n Analysis of coalesced communities"
Rscript analysis_community_coalescence.r 

echo -e "\nGENERATING REPORT"
cd writeup
#Run latex document with the write up
bash ~/CompileLaTeXstudio.sh Lechon_Pablo_CMEE_MSc_2020 1>../Sandbox/warnings.txt
cd ..
