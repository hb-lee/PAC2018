# PAC2018


### competition  

###### 1. Solver_pcg_mod.f90 - DJS

```fortran
do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      AX(i,j,bid) = A0 (i ,j ,bid)*X(i ,j ,bid) + &
                    AN (i ,j ,bid)*X(i ,j+1,bid) + &
                    AN (i ,j-1,bid)*X(i ,j-1,bid) + &
                    AE (i ,j ,bid)*X(i+1,j ,bid) + &
                    AE (i-1,j ,bid)*X(i-1,j ,bid) + &
                    ANE(i ,j ,bid)*X(i+1,j+1,bid) + &
                    ANE(i ,j-1,bid)*X(i+1,j-1,bid) + &
                    ANE(i-1,j ,bid)*X(i-1,j+1,bid) + &
                    ANE(i-1,j-1,bid)*X(i-1,j-1,bid)
   end do
   end do
```


git clone git@xxx  下载代码(只须在最开始下载，之后的操作不需要这个命令)   
git pull    拉取项目最新代码(同步远程仓到本地仓)     

本地修改代码后:   
git status 查看哪些修改文件;  
git add .     添加修改过的文件到本地仓库;   
git commit -m "xxx"    提交修改到本地仓;   
git push         提交本地仓的修改到远程仓库;   