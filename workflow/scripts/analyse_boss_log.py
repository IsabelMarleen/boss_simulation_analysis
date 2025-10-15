import re
import sys

def analyse_log(path, out_path):
    contr_line = 'control'
    boss_line = 'boss'
    all_batches_total = 0
    all_batches_reject = 0
    time_contr = 0
    time_boss = 0
    contr_dump = 1
    boss_dump = 1

    with open(out_path, "w") as o:
        o.write('cond,time,dump,total,base_total,unb,unb_ratio\n')
        with open(path, "r") as f:
            for line in f.readlines():
                time = re.findall(r'Next batch ---------------------------- # (\d+)', line)
                if len(time) != 0:
                    # contr_line = ','.join((contr_line, time[0], otu))
                    # boss_line = ','.join((boss_line, time[0], otu))
                    continue

                total = re.findall(r'got new batch of (\d+) reads', line)
                if len(total) != 0:
                    all_batches_total += int(total[0])
                    # contr_line = ','.join((contr_line, str(all_batches_total), str(0)))
                    # boss_line = ','.join((boss_line, str(all_batches_total), str(0)))
                    continue

                accept_rej = re.findall(r'.+ accepted (\d+), rejected (\d+)', line)
                if len(accept_rej) != 0:
                    all_batches_reject += int(accept_rej[0][1])
                    # contr_line = ','.join((contr_line, str(0), str(0)))+'\n'
                    # boss_line = ','.join((boss_line, str(all_batches_reject), str(all_batches_reject/all_batches_total)))+'\n'
                    continue

                tc = re.findall(r'.+ time control: (\d+)', line)
                if len(tc) != 0:
                    time_contr = float(tc[0])
                    continue

                tb = re.findall(r'.+ time boss-runs: (\d+)', line)
                if len(tb) != 0:
                    time_boss = float(tb[0])
                    continue

                dc = re.findall(r'.+ dump control #(\d+)', line)
                if len(dc) != 0:
                    contr_dump = int(dc[0])
                    contr_line = f"{contr_line},{time_contr},{contr_dump},{all_batches_total},{0},{0},{0}\n"
                    o.write(contr_line)
                    contr_line = 'control'
                    continue

                # o.write('cond,time,dump,total,base_total,unb,unb_ratio\n')
                db = re.findall(r'.+ dump boss #(\d+)', line)
                if len(db) != 0:
                    boss_dump = int(db[0])
                    boss_line = f"{boss_line},{time_boss},{boss_dump},{all_batches_total},{0},{all_batches_reject},{all_batches_reject/all_batches_total}\n"
                    o.write(boss_line)
                    boss_line = 'boss'
                    continue


if __name__ == "__main__":
    # input to this script is the logfile produced by the run
    analyse_log(path=sys.argv[1], out_path=sys.argv[2])