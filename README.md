# 下载
  ```
  git clone git@192.168.50.60:lanwangwei/dwave.git
  ```

# 安装
  ```
  cd dwave 
  conda-env create -n myenv -f=environment.yml
  cd ./dimod-main/
  python setup.py install
  cd ../dwave-samplers/
  python setup.py install
  ```

# 测试

  进入dwave文件夹
  ```
  python test.py
  ```
  output : 

>     0  1  2  3  4 energy num_oc.
> 0  0  1  1  0  1   -5.0       1
> 1  1  0  0  1  1   -5.0       1
> 2  1  0  0  1  0   -5.0       1
> 3  1  0  0  1  1   -5.0       1
> 4  1  0  0  1  1   -5.0       1
> 5  0  1  1  0  0   -5.0       1
> 6  0  1  1  0  1   -5.0       1
> 7  0  1  1  0  1   -5.0       1
> 8  1  0  0  1  0   -5.0       1
> 9  1  0  0  1  1   -5.0       1
> ['BINARY', 10 rows, 10 samples, 5 variables]
