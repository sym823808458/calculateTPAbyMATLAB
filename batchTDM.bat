for /f %%i in ('dir *.fchk /b') do (
Multiwfn %%i < Step_1.txt1
rename transdipmom.txt %%~ni.txt
)
//不想要输出过程，multiwfn自带输出文件的这样写