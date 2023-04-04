#!/bin/bash

#    # Variables
#    HOST=10.32.153.46
#    PORT=8889
#    USER=ovillegas
#    JUPYTER_PORT=8889
#    JUPYTER_PASSWORD="contrasena_de_jupyterlab"
#    
#    # Abrir el túnel SSH
#    ssh -N -L ${JUPYTER_PORT}:localhost:8888 ${USER}@${HOST} -p ${PORT} &
#    
#    # Iniciar Jupyter Lab
#    JUPYTER_TOKEN=`python -c "from notebook.auth import passwd; print(passwd('${JUPYTER_PASSWORD}'))"`
#    jupyter lab --ip=127.0.0.1 --port=${JUPYTER_PORT} --no-browser --NotebookApp.token=${JUPYTER_TOKEN}
#    
#    
#!/bin/bash
local=localhost:8889
remoto=8888
user=ovillegas
ip_servidor=10.32.153.46

ssh -L $local:localhost:$remoto $user@$ip_servidor -N &
ssh_pid=$!

#    ssh -N -L localhost:8889:localhost:8888 ovillegas@10.32.153.46
#    
#    ssh_pid=$!
#    
#    jupyter lab --no-browser --port=$puerto_remoto &
#    jupyter_pid=$!
#    
#    sleep 5 # Espera a que JupyterLab se inicie
#    
#    xdg-open "http://localhost:$puerto_local"
#    
#    read -p "Presione enter para detener el túnel SSH y JupyterLab" dummy
#    
#    kill $ssh_pid
#    kill $jupyter_pid
