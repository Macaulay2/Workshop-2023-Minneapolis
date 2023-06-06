## How to git? (simplified manual)
* Get a github account.
* The best practice (for github) is to set up an SSH-key: https://github.com/settings/keys
  (read https://help.github.com/articles/generating-an-ssh-key/ for instructions)
  * Clone a repo:
  **`git clone git@github.com:antonleykin/M2.git`**
  * Switch to your branch (if necessary; by default YOURBRANCH=_master_):
  `git checkout YOURBRANCH`
  * Change directory to the directory of interest (inside the cloned repo):
  `cd M2/M2/Macaulay2/packages`
  * Your "cycle" could be...
    * Pull in the updates from repo: **`git pull`**
      (do this often -- this updates your local repo with the changes others made; pull right before a commit -- that will save you a lot of trouble)
    * Commit your local edits: **`git commit -am "message for the log"`**
	  (this submits your updates to all files that are tracked by the repo; if you need the repo to track an additional file use `git add FILENAME` before committing)
    * Push the local commits: **`git push`** (this "sends" your changes to the remote repo)
    * At any time see the status (in particular, it tells you which branch you are on): **`git status`**
* Interactive tools:
    [magit](https://magit.vc/) (for emacs),
    GitHub Desktop,
    VSCode, ..., even Matlab.
* Conflict resolution: try to avoid, but...
  * the most typical conflict occurs when two people commit changes to the same part of a file: see [merge](https://help.github.com/articles/resolving-a-merge-conflict-from-the-command-line/) conflicts resolution.
  * **`git stash push`** and **`git stash pop`** could be used to "hide" and then reapply your _uncommitted_ changes: that helps in a scenario when `git pull` complains that you made changes (most likely to a different part) of the file that somebody else changed.

