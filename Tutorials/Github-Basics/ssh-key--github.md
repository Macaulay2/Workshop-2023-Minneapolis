# ssh key for GitHub (M2 version)

On your machine:
```
$ ssh-keygen -t ed25519 -C "your_email@example.com"
$ eval "$(ssh-agent -s)"
$ ssh-add ~/.ssh/id_ed25519
```
If you cloned you repo with `https` (as opposed to ssh), e.g. you see
```
$ git remote -v
origin  https://github.com/Macaulay2/Workshop-2023-Minneapolis.git (fetch)
origin  https://github.com/Macaulay2/Workshop-2023-Minneapolis.git (push)
```
switch over to `ssh` by
```
git remote set-url origin git@github.com:Macaulay2/Workshop-2023-Minneapolis.git
```

On github:
* [read this](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account#adding-a-new-ssh-key-to-your-account)
* follow this instructions
  