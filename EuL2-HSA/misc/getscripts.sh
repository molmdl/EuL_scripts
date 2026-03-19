#!/bin/bash
for i in in py sh txt gpi dat mdp inp ; do 
	find . -name "*.${i}*" -prune -type f >> scripts
done
echo ./scripts >> scripts

