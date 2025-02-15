# run_log_file.py
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import traceback

from .download import download_sra_to_variable, set_manager

def display_value(val):
    return str(val) if val is not None else ""

def convert_symbol(symbol):
    if symbol == "[green]✓[/green]":
        return "OK"
    elif symbol == "[red]✗[/red]":
        return "Failed"
    elif symbol == "[yellow]?[/yellow]":
        return "?"
    return symbol

def handle_future_exception(future, sra_id, log_file, updates):
    """
    This callback is fired immediately when the future completes, either successfully or with an exception.
    If there's an exception, we log it directly (without waiting for other processes to finish).
    """
    exc = future.exception()
    if exc is not None:
        updates[sra_id]["status"] = "Error"
        updates[sra_id]["downloaded_symbol"] = "[red]✗[/red]"
        with open(log_file, "a") as lf:
            lf.write(f"\nException in process for {sra_id}:\n")
            traceback.print_exc(file=lf)
        traceback.print_exc()

def run_with_log_file(sra_ids, threads, hmm_file, min_length, evalue, output_dir, fix_circles, log_file, save_all):
    manager = multiprocessing.Manager()
    set_manager(manager)  # Set the global manager for download.py
    updates = manager.dict()
    done_events = []

    executor = ProcessPoolExecutor(max_workers=threads)
    for sra_id in sra_ids:
        done_event = manager.Event()
        done_events.append(done_event)
        future = executor.submit(
            download_sra_to_variable,
            sra_id,
            updates,
            done_event,
            hmm_file,
            min_length,
            evalue,
            output_dir,
            fix_circles,
            save_all
        )
        future.add_done_callback(
            lambda fut, sid=sra_id: handle_future_exception(fut, sid, log_file, updates)
        )

    headers = ["name", "initial size", "status", "downloaded", "annotated", "HMMsearched", "contigs", "proteins", "hits"]
    with open(log_file, 'w') as lf:
        lf.write("\t".join(headers) + "\n")

        while True:
            completed = sum(1 for e in done_events if e.is_set())
            total = len(sra_ids)

            lf.seek(0)
            lf.truncate()
            lf.write("\t".join(headers) + "\n")

            updates_copy = list(updates.items())
            for sra_id, data in updates_copy:
                row = [
                    data["name"],
                    "",  # 'initial size' placeholder; adjust if needed
                    data["status"],
                    convert_symbol(data["downloaded_symbol"]),
                    convert_symbol(data["annotated_symbol"]),
                    convert_symbol(data["searched_symbol"]),
                    display_value(data.get("num_contigs")),
                    display_value(data.get("num_proteins")),
                    display_value(data.get("num_hits"))
                ]
                lf.write("\t".join(row) + "\n")

            lf.write(f"\nCompleted: {completed}/{total}\n")
            lf.flush()

            if completed == total:
                successful = []
                failed = []
                no_hits = []
                for sra_id, data in updates.items():
                    symbol = data.get("downloaded_symbol", "")
                    hits = data.get("num_hits", 0) or 0
                    if symbol == "[red]✗[/red]":
                        failed.append(sra_id)
                    elif symbol == "[green]✓[/green]":
                        if hits > 0:
                            successful.append(sra_id)
                        else:
                            no_hits.append(sra_id)
                    else:
                        no_hits.append(sra_id)

                lf.write("\nProcessing complete:\n")
                lf.write(f"Total SRAs processed: {len(sra_ids)}\n")
                lf.write(f"Successful (with hits): {len(successful)}\n")
                lf.write(f"No hits found: {len(no_hits)}\n")
                lf.write(f"Failed: {len(failed)}\n")
                lf.flush()
                break