
bash ./conduction/Postdamage_NewPlug.sh "11"
sleep 3600
bash ./conduction/Postdamage_sub_NewPlug.sh "11"
sleep 7200
bash ./conduction/Predamage_NewPlug.sh "11"
sleep 1800
bash ./conduction/FeymannKacs_prepare_NewPlug.sh "11"
sleep 120

bash ./conduction/FeymannKacs_simulate_NewPlug.sh "11"
