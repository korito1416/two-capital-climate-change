
bash ./conduction/Postdamage_NewPlug.sh "9"
sleep 3600
bash ./conduction/Postdamage_sub_NewPlug.sh "9"
sleep 3600
bash ./conduction/Predamage_NewPlug.sh "9"
sleep 1800
bash ./conduction/FeymannKacs_prepare_NewPlug.sh "9"
sleep 120

bash ./conduction/FeymannKacs_simulate_NewPlug.sh "9"

