#!/bin/bash


rsync -avz --progress --ignore-times xchi@seadragon:~/intern_backfilling/results/ ~/Desktop/Github_Local/backfilling/results/


echo "DONE Downloading Jobs."

