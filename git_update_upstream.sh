#!/bin/bash

## Update your local repo from SWAG repo. Creates upstream if necessary 

echo "Check if ustream exists"
git ls-remote --exit-code upstream
if test $? != 0; then
	echo "...NO"
	echo "Create upstream"
    git remote add upstream https://github.com/swag-kaust/TD.git
    echo "...OK"
    echo
else
	echo "...OK"
fi
echo

## Fetch updates from the remote repo
echo "Fetch"
git fetch upstream
if test $? = 0; then
echo "...OK"
else
	echo "ERROR"
fi
echo

## then: (like "git pull" which is fetch + merge)
echo "Merge master"
 git merge upstream/master master
 if test $? = 0; then
echo "...OK"
else
	echo "ERROR"
fi
## orreplay your local work on top of the fetched branch
## like a "git pull --rebase"
#git rebase upstream/master


echo
echo "DONE. Up to date"
