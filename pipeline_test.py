#!/home/chenzh/.conda/envs/cfRNA/bin/python3
#coding=utf-8
import argparse
import re
import os
import xml.etree.cElementTree as ET
import yaml
import subprocess
import queue
import time
import sys

# config = '/project/personal/chenzh/cfRNA/pipeline/config_test.yaml'
# xml = '/project/personal/chenzh/cfRNA/pipeline/pipeline_test.xml'

###临时默认变量：
#接头。默认为illumina
Primer_R1 = 'AGATCGGAAGAG'
Primer_R2 = 'AGATCGGAAGAG'

def parameter(): 
    with open(config,'r') as yaml_config:
        global conf_yaml,Project,Project_id,Group,Diff,Outdir,Sample,Anno_dir,Speices
        conf_yaml = yaml.load(yaml_config,Loader=yaml.Loader)    #Loader参数后面的yaml.Loader是载入模式。
        Project = conf_yaml['Project']
        Project_id = conf_yaml['Project_id']
        Speices = conf_yaml['Speices']
        Sample = conf_yaml['Sample']
        Group = conf_yaml['Group']
        Anno_dir = conf_yaml['Anno_dir']
        method = conf_yaml['method']
        Outdir = conf_yaml['Outdir']
        if 'Diff' in conf_yaml:
            Diff = conf_yaml['Diff']
        else:
            Diff = 'NA'
    # print(Project,Project_id,Speices,Group,Diff,Outdir,Anno_dir,sep='\t')  #查看打印出来的yaml.config输出


def write_shell(shell_file,content):
    with open(shell_file,'w') as shell_file_fp:
        content_final = re.sub("\n\s+","\n",content.lstrip())  ##这里是去除写入的开头的空格
        shell_file_fp.write(content_final)


def shell_del():
    tree = ET.parse(xml)   #解析xml文件
    root = tree.getroot()
    Task_dic={}  #任务字典，key为步骤序号，value为任务sh  
    for task in root.findall('Task'):  #Element.findall():仅查找带有标签的元素，这些元素是当前元素的直接子元素。
        shell_content = task.find('Coding').text
        shell_file_name = task.find('Coding_file').text
        step_id = re.search(r'shell_.*?_(\d+).*',shell_file_name).group(1).lstrip('0')
        if task.find('run_step').text == 'QC_MAP':   #判断执行哪些步骤以知道下面循环读取ymal文件的哪些内容
            for Sample_id,fq_path in Sample.items():
                data_R1 = fq_path[0]
                data_R2 = fq_path[1]
                content = shell_content.format(data_raw_R1=data_R1,data_raw_R2=data_R2,Project=Project,Sample_id=Sample_id,Outdir=Outdir,Anno_dir=Anno_dir,Speices=Speices,Primer_R1=Primer_R1,Primer_R2=Primer_R2)
                shell_name = shell_file_name.format(Project=Project,Sample_id=Sample_id,Outdir=Outdir)
                if not Task_dic.get(step_id):
                    Task_dic[step_id]=[]
                Task_dic[step_id].append(shell_name)
                write_shell(shell_name,content)
            try:
                shell_content = task.find('Coding2').text
                content = shell_content.format(Project=Project,Outdir=Outdir)
                shell_name = shell_file_name.format(Project=Project,Sample_id="All",Outdir=Outdir)
                step_id=f'{step_id}_2'
                if not Task_dic.get(step_id):
                    Task_dic[step_id]=[]                
                Task_dic[step_id].append(shell_name)
                write_shell(shell_name,content)
            except:
                pass

        if task.find('run_step').text == 'Count': 
            content = shell_content.format(Project=Project,Outdir=Outdir,Anno_dir=Anno_dir,Speices=Speices)
            shell_name = shell_file_name.format(Project=Project,Outdir=Outdir)
            if not Task_dic.get(step_id):
                Task_dic[step_id]=[]
            Task_dic[step_id].append(shell_name)
            write_shell(shell_name,content)

        if task.find('run_step').text == 'Diff':
            # for sample_num in Group.values():
            #     if len(sample_num)>1:
            #         Judge="DEseq2"
            #     else:
            #         Judge="edgeR"
            for DEG_group in Diff:
                Treat=DEG_group.split("_vs_")[0]
                Control=DEG_group.split("_vs_")[1]
                if len(Group[Treat]) > 1 and len(Group[Control])>1:
                    Judge="DEseq2"
                else:
                    Judge="edgeR"
                Treat_list=','.join(Group[Treat])
                Control_list=','.join(Group[Control])
                content=shell_content.format(Project=Project,Outdir=Outdir,Treat=Treat,Treat_list=Treat_list,Control=Control,Control_list=Control_list,Judge=Judge,Diff=DEG_group)
                shell_name=shell_file_name.format(Project=Project,Outdir=Outdir,Diff=DEG_group)
                if not Task_dic.get(step_id):
                    Task_dic[step_id]=[]
                Task_dic[step_id].append(shell_name)
                write_shell(shell_name,content)     
    return(Task_dic)

        




if __name__ == '__main__':
    doc='''
    ##########
    chenzh: cfRNA测试
    time: 2023.05.09
    version: 1.0
    ##########
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description='The pipeline for cfRNA',prog='cfRNA-seq program',epilog=doc)
    parser.add_argument('--cfg',help='Config file of Data path and Custom Infomation')
    parser.add_argument('--job',help='Number of threads running the program(default %(default)s)',default=6,type=int)
    parser.add_argument('--xml',help='ChIPseq.xml (default %(default)s)',default='/home/chenzh/pipeline/cfRNA/pipeline_test.xml')
    parser.add_argument('--str',help='str program (default %(default)s)',default=0)
    parser.add_argument('--end',help='end program (default %(default)s)',default=20)
    parser.add_argument('--inexe',nargs='?',default=1,help='0 表示检查输出脚本，1表示生成脚本并运行流程，也可以不写')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)

    config=args.cfg
    jobs_num=int(args.job)
    xml=args.xml
    str=int(args.str)
    end=int(args.end)
    inexe=int(args.inexe)

    parameter()

    #创建文件夹
    mkdir_dir = os.path.join('%s/%s/{01.Data,02.QC,03.MAP,04.EXP,05.Diff}'%(Outdir,Project))
    mkdir_dir2 = os.path.join('%s/Task/log'%Outdir)
    cmd_1 = 'mkdir -p %s' % (mkdir_dir)
    cmd_2 = 'mkdir -p %s' % (mkdir_dir2)
    os.system(cmd_1)
    os.system(cmd_2)


    Task_dic=shell_del()

    if inexe:
        running_tasks=[]
        max_jobs=jobs_num
        #投递任务函数
        # def submit_jobs(job_dic,max_jobs):
        task_que=queue.Queue()

        #将任务添加进队列
        for key,value in Task_dic.items():
            # print(key,":",value)
            task_que.put((key,value))
            # print(task_que.get())
        
        while not task_que.empty() or running_tasks:
            task_id,task_list=task_que.get()
            task_id=task_id.split("_")[0]
            if int(task_id) in range(str,end):
                while len(task_list)!=0:
                    for task in task_list:
                        if len(running_tasks) < max_jobs and (task_id,task) not in running_tasks:
                            cmd=('sh %s &'%task)
                            log_name=task.split('/')[-1]
                            with open(f'{Outdir}/Task/log/{log_name}.o','w') as o:
                                with open(f'{Outdir}/Task/log/{log_name}.e','w') as e:
                                    subprocess.run(cmd,shell=True,stdout=o,stderr=e)
                            running_tasks.append((task_id,task))
                    # print(running_tasks)
                        else:
                            break
                        
                    #检查正在运行的任务是否小于并发数，且队列中还有任务，则将任务从队列中取出并运行：
                    for rm_task in reversed(running_tasks):
                        rm_task_id,rm_task_data = rm_task
                        output=subprocess.run('ps -fx |grep %s |wc -l'%rm_task_data,shell=True,stdout=subprocess.PIPE).stdout.decode('utf-8')
                        if not '3\n' in output:
                            print(f'Task {rm_task_data} finished')
                            running_tasks.remove(rm_task)
                            task_list.remove(rm_task_data)
                    
                    time.sleep(10)
