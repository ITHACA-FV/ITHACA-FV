**Developer Instructions**

In order to start developing the code and to start adding your contributions to the repository please follow these steps. If it is the first time you start adding contributions to the ITHACA-FV library follow these steps. From point number 1) If you have already contributed to the project you can start directly from point 4). However, before creating the development branch that contains the new feature you want to add, remember to pull the last changes from the mathlab remote master branch:
```
git pull mathlab master
```
1) Make a fork of the public repository onto your gihub account,  https://github.com/mathlab/ITHACA-FV and let's suppose from now on that your user account is **pincopallino**. 
2) Once you have forked the repository you can clone it to your pc:
```
git clone https://github.com/pincopallino/ITHACA-FV
```
3) Go the cloned folder and add the remote URL of the mathlab version of the ITHACA-FV repository
```
cd ITHACA-FV
git remote add mathlab https://github.com/mathlab/ITHACA-FV
```
4) Create a development branch into your local repository
```
git checkout -b development
```
5) Start implementing changes or adding new functionalities
6) Add changed files or new files for the next commit
```
git add NEWFILE1.H
git add NEWFILE2.C
```
7) Prepare a commit 
```
git commit -m "meaningfull commit message"
```
8) Push the development branch commit to your remote repository
```
git push origin development
```
9) Perform a pull request from the github webpage of ITHACA-FV, the mathlab one at https://github.com/mathlab/ITHACA-FV. You will have to click " New pull request" and then on "compare across forks". You will have to compare the master branch of the mathlab version of ITHACA "mathlab/ITHACA-FV/master" with the development branch of your version of ITHACA-FV "pincopallino/ITHACA-FV/development". Once you have performed the pull request you will have to wait the response of the ITHACA-FV administrators. 
10) In case there is something to change you will have to implement corrections asked by the administrators and perform a new push. In case there are no corrections to be done you don't have to consider this point.
```
git add "modifiedfile.H"
git commit -m "meaningfull commit message"
git push origin development
```
11) Once the pull request has been accepted you can checkout on the master branch of your local repo, pull the last changes from the mathlab remote master branch, delete the development branch locally and remotely and finally push the local changes onto your remote repository.
```
git checkout master
git pull mathlab master
git branch -d development
git push origin --delete development
git push origin master
```

**HAVE FUN CODING!!!**

