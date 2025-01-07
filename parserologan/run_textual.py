import sys
import asyncio
import threading
from io import StringIO
from concurrent.futures import ThreadPoolExecutor
from textual.app import App, ComposeResult
from textual.widgets import DataTable, Header, Footer, Static, RichLog
from textual.containers import Vertical
from textual.reactive import reactive
from .download import download_sra_to_variable

def display_value(val):
    return str(val) if val is not None else ""

def run_with_textual_interface(sra_ids, threads, hmm_file, min_length, evalue, output_dir, fix_circles,save_all=False):
    class LoganParsero(App):
        CSS = """
        Screen {
            layout: vertical;
        }
        """

        sra_ids = []
        updates = {}
        done_events = []
        sra_to_row = {}
        done_flag = threading.Event()
        executor = ThreadPoolExecutor()
        completed = reactive(0)

        def __init__(self, sra_ids, threads, hmm_file, min_length, evalue, output_dir, fix_circles=False):
            super().__init__()
            self.sra_ids = sra_ids
            self.threads = threads
            self.hmm_file = hmm_file
            self.min_length = min_length
            self.evalue = evalue
            self.output_dir = output_dir
            self.fix_circles = fix_circles
            self.save_all = save_all
            self._stdout_buffer = StringIO()

        def compose(self) -> ComposeResult:
            yield Header()
            yield Vertical(
                Static("SRA Download, Annotation, and HMM Search Progress", id="title"),
                DataTable(id="table"),
                Static("", id="status"),
                RichLog(id="log", highlight=True)
            )
            yield Footer()

        async def on_mount(self):
            sys.stdout = self._stdout_buffer

            table = self.query_one("#table", DataTable)
            self.name_col = table.add_column("name", width=25)
            self.size_col = table.add_column("initial size")
            self.status_col = table.add_column("status")
            self.downloaded_col = table.add_column("downloaded")
            self.annotated_col = table.add_column("annotated")
            self.searched_col = table.add_column("HMMsearched")
            self.contig_col = table.add_column("contigs")
            self.protein_col = table.add_column("proteins")
            self.hits_col = table.add_column("hits")

            for sra_id in self.sra_ids:
                done_event = threading.Event()
                self.done_events.append(done_event)
                self.executor.submit(
                    download_sra_to_variable,
                    sra_id,
                    self.updates,
                    done_event,
                    self.hmm_file,
                    self.min_length,
                    self.evalue,
                    self.output_dir,
                    self.fix_circles,
                    self.save_all
                )

            self.set_interval(0.5, self.refresh_table)

        def refresh_table(self):
            table = self.query_one("#table", DataTable)
            completed = sum(1 for e in self.done_events if e.is_set())
            total = len(self.sra_ids)

            status = self.query_one("#status", Static)
            status.update(f"Completed: {completed}/{total}")

            updates_copy = list(self.updates.items())

            for sra_id, data in updates_copy:
                if sra_id not in self.sra_to_row:
                    row_key = table.add_row(
                        data["name"],
                        data["total_size"],
                        data["status"],
                        data["downloaded_symbol"],
                        data["annotated_symbol"],
                        data["searched_symbol"],
                        display_value(data.get("num_contigs")),
                        display_value(data.get("num_proteins")),
                        display_value(data.get("num_hits"))
                    )
                    self.sra_to_row[sra_id] = row_key
                else:
                    row_key = self.sra_to_row[sra_id]
                    table.update_cell(row_key, self.name_col, data["name"])
                    table.update_cell(row_key, self.size_col, data["total_size"])
                    table.update_cell(row_key, self.status_col, data["status"])
                    table.update_cell(row_key, self.downloaded_col, data["downloaded_symbol"])
                    table.update_cell(row_key, self.annotated_col, data["annotated_symbol"])
                    table.update_cell(row_key, self.searched_col, data["searched_symbol"])
                    table.update_cell(row_key, self.contig_col, display_value(data.get("num_contigs")))
                    table.update_cell(row_key, self.protein_col, display_value(data.get("num_proteins")))
                    table.update_cell(row_key, self.hits_col, display_value(data.get("num_hits")))

            log = self.query_one("#log", RichLog)
            output = self._stdout_buffer.getvalue()
            if output:
                self._stdout_buffer.seek(0)
                self._stdout_buffer.truncate(0)
                for line in output.splitlines():
                    log.write(line)

            if completed == total and not self.done_flag.is_set():
                self.done_flag.set()
                asyncio.create_task(self.finish_up())

        async def finish_up(self):
            await asyncio.sleep(1)
            successful = []
            failed = []
            no_hits = []
            for sra_id, data in self.updates.items():
                symbol = data.get("downloaded_symbol", "")
                hits = data.get("num_hits", 0) or 0
                if "[red]✗[/red]" in symbol:
                    failed.append(sra_id)
                elif "[green]✓[/green]" in symbol:
                    if hits > 0:
                        successful.append(sra_id)
                    else:
                        no_hits.append(sra_id)
                else:
                    no_hits.append(sra_id)

            print("\nProcessing complete:")
            print(f"Total SRAs processed: {len(self.sra_ids)}")
            print(f"Successful (with hits): {len(successful)}")
            print(f"No hits found: {len(no_hits)}")
            print(f"Failed: {len(failed)}")

            await self.action_quit()

    LoganParsero(sra_ids, threads, hmm_file, min_length, evalue, output_dir, fix_circles).run()
