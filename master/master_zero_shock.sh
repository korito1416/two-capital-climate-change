#! /bin/bash


## Table 1 and Figure 5
bash ./conduction/Postdamage.sh "true" "false" "false" "false" "false"
sleep 1200
bash ./conduction/Postdamage_sub.sh "true" "false" "false" "false" "false"
sleep 1200
bash ./conduction/Predamage.sh "true" "false" "false" "false" "false"
sleep 300
bash ./conduction/ZeroShockTrajectories_simulate.sh "true" "false" "false" "false" 
sleep 120
bash ./conduction/ZeroShockTrajectories_plot.sh "true" "false" "false" "false" 

## Figure 4
bash ./conduction/Postdamage.sh "false" "false" "false" "true" "false"
sleep 3000
bash ./conduction/Postdamage_sub.sh "false" "false" "false" "true" "false"
sleep 3000
bash ./conduction/Predamage.sh "false" "false" "false" "true" "false"
sleep 300
bash ./conduction/ZeroShockTrajectories_simulate.sh "false" "false" "false" "true" 
sleep 120
bash ./conduction/ZeroShockTrajectories_plot.sh "false" "false" "false" "true" 


