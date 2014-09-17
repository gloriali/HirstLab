Script to retrieve metadata from LIMS
=====================================
* pwd is gin pwd.     
* API script

```
export PYTHONPATH=PYTHONPATH:/home/sgakkhar/work/admin/mypylibs
export asGuaRMyqydTHdqgLEqYaVdyOFBVKugagdPJWCQaq="GIN_PWD"
/home/sgakkhar/Documents/mypylibs/limsRefactored.py <libID> --simple
```

* for indexed library         

```
/home/sgakkhar/Documents/mypylibs/limsRefactored.py <libID> --rearray-target
```       

* To force it to promt me for a password

```
unset asGuaRMyqydTHdqgLEqYaVdyOFBVKugagdPJWCQaq 
````

* Check username that script will use         

```
echo $USER 
```      
