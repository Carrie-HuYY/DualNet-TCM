import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import threading
import os
import main
from datetime import datetime


def run_analysis():
    try:
        # 获取所有输入参数
        search_type = search_type_var.get()
        search_name = search_name_entry.get().strip()  # 获取搜索名称并去除首尾空格
        disease_name = disease_name_entry.get().strip()

        # 验证必要参数
        if not search_name:
            messagebox.showerror("错误", "搜索名称不能为空!")
            return
        if not disease_name:
            messagebox.showerror("错误", "疾病名称不能为空!")
            return

        try:
            score = int(score_entry.get())
            target_max_number = int(target_max_number_entry.get())
            report_number = int(report_number_entry.get())
            interaction_number = int(interaction_number_entry.get())
        except ValueError:
            messagebox.showerror("输入错误", "请输入有效的数值!")
            return

        out_graph = out_graph_var.get()
        out_for_cytoscape = out_for_cytoscape_var.get()
        out_for_excel = out_for_excel_var.get()
        research_status_test = research_status_test_var.get()
        safety_research = safety_research_var.get()
        path = output_path_entry.get().strip()

        # 确保输出路径存在
        if not path:
            path = "results"
        os.makedirs(path, exist_ok=True)

        # 添加时间戳到日志文件名
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = os.path.join(path, f"analysis_log_{timestamp}.txt")

        # 禁用按钮防止重复点击
        start_button.config(state=tk.DISABLED)
        status_label.config(text="分析进行中...", fg="blue")

        # 在新线程中运行分析函数
        def analysis_thread():
            try:
                with open(log_file, "w") as f:
                    f.write(f"开始分析于 {timestamp}\n")
                    f.write(f"参数: SearchType={search_type}, SearchName={search_name}, DiseaseName={disease_name}\n")

                # 将搜索名称转换为列表
                search_name_list = [search_name]  # 将字符串转换为列表

                # 调用主分析函数
                main.TCM_VOTER(
                    SearchType=search_type,
                    SearchName=search_name_list,  # 传递列表
                    DiseaseName=disease_name,
                    target_max_number=target_max_number,
                    report_number=report_number,
                    interaction_number=interaction_number,
                    score=score,
                    out_graph=out_graph,
                    out_for_cytoscape=out_for_cytoscape,
                    out_for_excel=out_for_excel,
                    research_status_test=research_status_test,
                    safety_research=safety_research,
                    path=path
                )

                # 检查生成的文件
                generated_files = []
                if out_graph and os.path.exists(os.path.join(path, "visualization.html")):
                    generated_files.append("visualization.html")
                if out_for_cytoscape and os.path.exists(os.path.join(path, "for_cytoscape.txt")):
                    generated_files.append("for_cytoscape.txt")
                if out_for_excel and os.path.exists(os.path.join(path, "results.xlsx")):
                    generated_files.append("results.xlsx")

                if generated_files:
                    message = f"分析完成！生成的文件: {', '.join(generated_files)}\n路径: {os.path.abspath(path)}"
                    status_label.config(text="分析完成", fg="green")
                else:
                    message = "分析完成，但未生成任何输出文件。请检查参数设置。"
                    status_label.config(text="无输出文件", fg="orange")

                with open(log_file, "a") as f:
                    f.write(f"{message}\n")

                messagebox.showinfo("完成", message)

            except Exception as e:
                error_msg = f"分析过程中发生错误：{str(e)}"
                status_label.config(text="分析失败", fg="red")
                with open(log_file, "a") as f:
                    f.write(f"错误: {error_msg}\n")
                messagebox.showerror("错误", error_msg)
            finally:
                start_button.config(state=tk.NORMAL)

        threading.Thread(target=analysis_thread, daemon=True).start()

    except Exception as e:
        messagebox.showerror("错误", f"发生意外错误: {str(e)}")
        start_button.config(state=tk.NORMAL)


def select_output_path():
    path = filedialog.askdirectory()
    if path:  # 用户没有取消选择
        output_path_entry.delete(0, tk.END)
        output_path_entry.insert(0, path)


# 创建主窗口
root = tk.Tk()
root.title("TCM_VOTER 网络药理学分析工具")

# 设置网格布局
root.columnconfigure(1, weight=1)
for i in range(14):
    root.rowconfigure(i, pad=5)

# 搜索类型 - 改为数字选项
tk.Label(root, text="搜索类型(0-4):").grid(row=0, column=0, sticky="e")
search_type_var = tk.IntVar()
search_type_combobox = ttk.Combobox(root, textvariable=search_type_var,
                                    values=[0, 1, 2, 3, 4])
search_type_combobox.grid(row=0, column=1, sticky="ew", padx=5)
search_type_combobox.current(0)  # 设置默认选项

# 添加搜索类型说明标签
type_explanation = tk.Label(root, text="0:证候 1:方剂 2:中药 3:成分 4:靶点", fg="gray")
type_explanation.grid(row=0, column=2, sticky="w", padx=5)

# 搜索名称
tk.Label(root, text="搜索名称:").grid(row=1, column=0, sticky="e")
search_name_entry = tk.Entry(root)
search_name_entry.grid(row=1, column=1, sticky="ew", padx=5)

# 疾病名称
tk.Label(root, text="疾病名称:").grid(row=2, column=0, sticky="e")
disease_name_entry = tk.Entry(root)
disease_name_entry.insert(0, "cough")
disease_name_entry.grid(row=2, column=1, sticky="ew", padx=5)

# 分数阈值
tk.Label(root, text="分数阈值:").grid(row=3, column=0, sticky="e")
score_entry = tk.Entry(root)
score_entry.insert(0, "990")
score_entry.grid(row=3, column=1, sticky="ew", padx=5)

# 最大靶点数量
tk.Label(root, text="最大靶点数量:").grid(row=4, column=0, sticky="e")
target_max_number_entry = tk.Entry(root)
target_max_number_entry.insert(0, "70")
target_max_number_entry.grid(row=4, column=1, sticky="ew", padx=5)

# 报告数量
tk.Label(root, text="报告数量:").grid(row=5, column=0, sticky="e")
report_number_entry = tk.Entry(root)
report_number_entry.insert(0, "0")
report_number_entry.grid(row=5, column=1, sticky="ew", padx=5)

# 相互作用数量
tk.Label(root, text="相互作用数量:").grid(row=6, column=0, sticky="e")
interaction_number_entry = tk.Entry(root)
interaction_number_entry.insert(0, "0")
interaction_number_entry.grid(row=6, column=1, sticky="ew", padx=5)

# 输出路径
tk.Label(root, text="输出路径:").grid(row=7, column=0, sticky="e")
output_path_entry = tk.Entry(root)
output_path_entry.insert(0, "results")
output_path_entry.grid(row=7, column=1, sticky="ew", padx=5)
tk.Button(root, text="浏览", command=select_output_path).grid(row=7, column=2, padx=5)

# 输出选项
tk.Label(root, text="输出选项:").grid(row=8, column=0, sticky="ne")
output_frame = tk.Frame(root)
output_frame.grid(row=8, column=1, columnspan=2, sticky="w")

out_graph_var = tk.BooleanVar(value=True)
tk.Checkbutton(output_frame, text="生成图形", variable=out_graph_var).pack(anchor="w")
out_for_cytoscape_var = tk.BooleanVar(value=True)
tk.Checkbutton(output_frame, text="Cytoscape文件", variable=out_for_cytoscape_var).pack(anchor="w")
out_for_excel_var = tk.BooleanVar(value=True)
tk.Checkbutton(output_frame, text="Excel文件", variable=out_for_excel_var).pack(anchor="w")

# 研究选项
research_frame = tk.Frame(root)
research_frame.grid(row=11, column=1, columnspan=2, sticky="w")

research_status_test_var = tk.BooleanVar(value=True)
tk.Checkbutton(research_frame, text="研究状态测试", variable=research_status_test_var).pack(anchor="w")
safety_research_var = tk.BooleanVar(value=True)
tk.Checkbutton(research_frame, text="安全性研究", variable=safety_research_var).pack(anchor="w")

# 状态标签
status_label = tk.Label(root, text="准备就绪", fg="gray")
status_label.grid(row=12, column=1, sticky="w")

# 开始按钮
start_button = tk.Button(root, text="开始分析", command=run_analysis)
start_button.grid(row=13, column=1, pady=10, sticky="e")

# 运行主循环
root.mainloop()