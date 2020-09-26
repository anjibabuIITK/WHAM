#
# This script is a part of TASS Reweighting Pachage
# Commands to use git 
# 
# Authour: Anji Babu Kapakayala
#
# USAGE: chmod 777 GIT
#        ./GIT KEYWORD
#
#!/bin/bash
case "$1" in

  "pull"|"pl")
	git pull origin master;;
  "push"|"ps")
	git push origin master;;
  "status"|"S"|"s")
	git status ;;
  "add"|"A"|"a")
	git add *
	git status;;
  "commit"|"c"|"C")
      read -p "Write your message:" msg
	git commit -m "$msg"
	git status;;
  "clone")
      read -p "Write your message:" link
        git clone "$link"
	git status;;
  "rm")
      read -p "Enter Filename Carefully:" File
	git rm $File
	git status;;
  "*")
	echo "Give appropriate git command";;
esac
