#!/usr/bin/env zsh

if [[ ! -f .git/config ]]
then
	echo "must be run from the root dir of the repository"
	exit 1
fi

if [[ ! -d notebook ]]
then
	echo "creating notebook dir ..."
	mkdir notebook
fi

pushd notebook
	echo "creating .gitattributes in notebook dir"
	echo "*.ipynb filter=strip-notebook-output" > .gitattributes
popd

echo "creating .gitconfig rule to convert notebook"
cat << EOT > .gitconfig
[filter "strip-notebook-output"]
	clean = "jupyter nbconvert --clear-output --to=notebook --stdin --stdout --log-level=INFO"
EOT

echo "adding .gitconfig to the local repository"
git config --local include.path ../.gitconfig
