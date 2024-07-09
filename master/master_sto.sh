## Figure 8

bash ./conduction/StocasticTrajectories_simulate.sh "true" "false" 
sleep 1800


bash ./conduction/Postdamage.sh "false" "false" "false" "false" "true"
sleep 1800
bash ./conduction/Postdamage_sub.sh "false" "false" "false" "false" "true"
sleep 1800
bash ./conduction/Predamage.sh "false" "false" "false" "false" "true"
sleep 300
bash ./conduction/StocasticTrajectories_simulate.sh "false" "true" 

