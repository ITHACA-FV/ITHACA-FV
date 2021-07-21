## Developer Instructions

In order to start developing the code and to start adding your contributions to
the repository please follow these steps. In particular:
* if it is the first time you start adding contributions to the ITHACA-FV library
follow these steps from point number [1](#1);
* if you have already contributed to the project you can start directly from point number [4](#4).

However, before creating the development branch that contains the new feature
you want to add, remember to pull the last changes (see [here](https://help.github.com/articles/syncing-a-fork/) for more information).

If you intend to add some new classes or methods please follow the naming convention reported [here](https://gist.github.com/lefticus/10191322) and in particular use the CamelCase convention.

1. <div id="1">Create a new fork of the public repository
   [ITHACA-FV](https://github.com/mathlab/ITHACA-FV)
   ([here](https://help.github.com/articles/fork-a-repo) you can find an
    introduction about the fork procedure). 

2. Let's assume your Github username is YOUR-USERNAME, you can access to your
   own fork at the URL
   [github.com/YOUR-USERNAME/ITHACA-FV](https://github.com/YOUR-USERNAME/ITHACA-FV)
   and clone it using the favourite protocol (HTTPS is the easiest) using the
   command `git clone <url>`. In this example by digiting on the terminal:
   ```
   git clone "https://github.com/YOUR-USERNAME/ITHACA-FV"
   ```
   a new folder called ITHACA-FV is created.

3. Go the cloned folder and add the remote URL of the original ITHACA-FV
   repository (**not** the fork) using the following commands:
   ```
   cd ITHACA-FV
   git remote add mathlab "https://github.com/mathlab/ITHACA-FV"
   ```

4. <div id="4">Create a branch into your local repository. The new branch name
   should describe the new feature you have planned to insert in the
   repository. The command is:
   ```
   git checkout -b new_branch_name
   ```

5. Start implementing changes or adding new functionalities. Please ensure to
   write well-documented and properly formatted code. You can use the
   [`code_formatter.sh`](https://github.com/mathLab/ITHACA-FV/blob/master/code_formatter.sh)
   script (already in the repository) to automatically format the code
   according to the ITHACA style.

6. Add changed files or new files for the next commit by digiting:
   ```
   git add NEWFILE1.H
   git add NEWFILE2.C
   ```

7. Create the new commit. The commit message has to synthetically describe the
   new changes (see [here](https://chris.beams.io/posts/git-commit/) for an
   exhaustive discussion about commit messages).
   ```
   git commit -m "meaningfull commit message"
   ```

8. Before pushing your local changes in the Github fork, don't forget to sync
   your fork with the original repository. With this further step, you can
   avoid conflicts and maintain a clean history of the repository. To do this,
   you have to fetch the new changes, pull the `master` branch then adjust the
   history of your `new_branch_name` branch. Shortly:
   ```
   git checkout master
   git fetch mathlab
   git pull mathlab master
   git checkout new_branch_name
   git rebase master new_branch_name
   ```
   For a detailed discussion about rebasing and synching, you can read
   [this](https://help.github.com/articles/syncing-a-fork/).

9. Push the `new_branch_name` branch commit to your remote fork:
   ```
   git push -f origin new_branch_name
   ```

10. Open a new pull request (in [this guide](https://help.github.com/articles/creating-a-pull-request/)
    the basic steps to open a pull request are introduced). Please use a
    meaningful name and a synthetic description of the new implemented feature.
    Once you have performed the pull request you will have to wait the response
    of the ITHACA-FV administrators. 

11. In case there is something to change you will have to implement corrections
    asked by the administrators and perform a new push. In case there are no
    corrections to be done you don't have to consider this point.
    ```
    git add "modifiedfile.H"
    git commit -m "another meaningfull commit message"
    git push origin new_branch_name
    ```

12. Once the pull request has been accepted you have to synchronize your fork
    with the original repository (at this moment, only the original repository, 
    and not your fork, contains within the `master` branch the feature you
    implemented).
    As before digit (*NOTE*: the following commands delete the `new_branch_name`):
    ```
    git checkout master
    git fetch mathlab
    git pull mathlab master
    git branch -D new_branch_name
    git push origin --delete new_branch_name
    git push origin master
    ```
**PULL REQUESTS CONVENTIONS*

Use always tags to identify your pull request in order to make the life of maintainers easier:

See how a minor change to your commit message style can make you a better programmer.

Format: `<type>(<scope>): <subject>`

`<scope>` is optional

## Example

```
feat: add hat wobble
^--^  ^------------^
|     |
|     +-> Summary in present tense.
|
+-------> Type: chore, docs, feat, fix, refactor, style, or test.
```

## Type of Pull Requests

- `FEAT`: (new feature for the user, not a new feature for build script)
- `ENH`: (Enhance existing code without the inclusions of additional features)
- `FIX`: (bug fix for the user, not a fix to a build script)
- `DOCS`: (changes to the documentation)
- `STYLE`: (formatting, missing semi colons, etc; no production code change)
- `REFACTOR`: (refactoring production code, eg. renaming a variable)
- `TEST`: (adding missing tests, refactoring tests; no production code change)

**HAVE FUN CODING!!!**

