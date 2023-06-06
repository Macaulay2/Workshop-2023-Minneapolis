M-x shell

## Working alone

```
git clone git@github.com:Macaulay2/Workshop-2023-Minneapolis.git B
cd B
git checkout git-tutorial
```

## The basic "cycle"

To update (often!)
```
git pull
```
After making changes
```
git commit -a -m "I did it!"
git push
```

If you created FILENAME
```
git add FILENAME
```
Use `git status` often.

## Working with friends

Strategy 1. Edit separate files.

Strategy 2. Try to edit separate parts of the file, if you have to edit the same one.

Strategy 3. Be ready to resolve conflicts.

Scenario: something changes (on remote repo) in the file that I edit.
```
git stash push
git pull
git stash pop
```
Scenario: some change conflicts with my edits.
```
git status
```
Edit the files (e.g., FILENAME) with conflicts.
```
git add FILENAME
git commit -m "resolved conflict: yay!"
```

## It is a matter of style!

GitHub has many other tools: e.g., look at "Projects" and "Wiki".

How you use `git` depends on what works for you.


