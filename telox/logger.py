# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

import sys
import logging
import subprocess
import time
import resource
import shlex
import shutil
from pathlib import Path

class SmartFormatter(logging.Formatter):
    def format(self, record):
        original_fmt = self._style._fmt
        if getattr(record, 'bare', False):
            self._style._fmt = "%(message)s"
        elif hasattr(record, 'step_tag'):
            self._style._fmt = "[%(step_tag)s] %(message)s"
        result = super().format(record)
        self._style._fmt = original_fmt
        return result

class WorkflowLogger:
    def __init__(self, outdir: Path, prefix: str, debug: bool = False, mini_log: bool = False, total_steps: int = 1, start_step: int = 1, version: str = "0.0.0"):
        self.outdir = outdir
        self.prefix = prefix
        self.debug = debug
        self.mini_log = mini_log
        self.total_steps = total_steps
        self.version = version

        self.current_step_num = start_step - 1
        self.step_start_time = None
        self.start_global = time.perf_counter()

        log_file = outdir / f"{prefix}.log"

        self.logger = logging.getLogger()
        self.logger.setLevel(logging.DEBUG if debug else logging.INFO)
        self.logger.handlers = []

        file_fmt = SmartFormatter('[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        
        fh = logging.FileHandler(str(log_file), mode='w')
        fh.setFormatter(file_fmt)
        fh.setLevel(logging.DEBUG if debug else logging.INFO)
        self.logger.addHandler(fh)

        console_fmt = SmartFormatter('[%(levelname)s] %(message)s')
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(console_fmt)
        ch.setLevel(logging.DEBUG if debug else logging.INFO)
        self.logger.addHandler(ch)

    def _validate_command(self, command):
        if isinstance(command, str):
            executable = command.split()[0]
        else:
            executable = command[0]
        
        if not shutil.which(executable):
            raise ValueError(f"Required command '{executable}' not found in PATH.")
        return executable

    def log_header(self, cmd_args, input_file=None, ref_file=None, preset_name=None, motif=None):
        clean_args = ["telox"] + list(cmd_args)[1:]
        cmd_str = " ".join(shlex.quote(arg) for arg in clean_args)
        
        if self.mini_log:
            self.logger.info(f"\nTeloXplorer v{self.version}", extra={'bare': True})
            self.logger.info(f"Command: {cmd_str}", extra={'bare': True})
            return

        preset_display = f"{preset_name} ({motif})" if (preset_name and motif) else preset_name
        
        header_lines = [
            f"\nTeloXplorer v{self.version}",
            '-' * 50,
            f"* CMD    : {cmd_str}"
        ]

        rows = [
            ("Input", input_file),
            ("Ref", ref_file),
            ("Output", f"{self.outdir}/{self.prefix}"),
            ("Preset", preset_display)
        ]

        for label, value in rows:
            if value:
                header_lines.append(f"* {label:<6} : {value}")

        self.logger.info("\n".join(header_lines), extra={'bare': True})
        self.logger.info("",extra={'bare': True})

    def start_step(self, name: str):
        self.current_step_num += 1
        self.step_start_time = time.perf_counter()
        display_num = min(self.current_step_num, self.total_steps)

        if self.mini_log:
            try:
                self.logger.handlers[1].terminator = "" 
                self.logger.info(f"[{display_num}/{self.total_steps}] {name}...", extra={'bare': True})
            finally:
                self.logger.handlers[1].terminator = "\n"
        else:
            step_label = f"STEP {display_num}/{self.total_steps}"
            self.logger.info(f"{name}...", extra={'step_tag': step_label})

    def end_step(self):
        if self.step_start_time is None:
            self.logger.warning("end_step called before start_step!")
            return

        elapsed = time.perf_counter() - self.step_start_time
        self.step_start_time = None

        if self.mini_log:
            self.logger.info(f" done. ({elapsed:.1f}s)", extra={'bare': True})
        else:
            self.logger.info(f"Done! ({elapsed:.1f}s)\n")

    def skip_step(self, name: str):
        self.current_step_num += 1
        display_num = min(self.current_step_num, self.total_steps)

        msg = "Skipping"

        if self.mini_log:
            self.logger.info(f"[{display_num}/{self.total_steps}] {msg}", extra={'bare': True})
        else:
            step_label = f"STEP {display_num}/{self.total_steps}"
            self.logger.info(msg, extra={'step_tag': step_label})

    def info(self, msg, **kwargs):
        if self.mini_log:
            fh = self.logger.handlers[0]
            record = logging.LogRecord(
                name=self.logger.name, 
                level=logging.INFO, 
                pathname="", 
                lineno=0, 
                msg=msg, 
                args=None, 
                exc_info=None
            )
            fh.emit(record)
        else:
            self.logger.info(msg, **kwargs)


    def success(self, module="teloxplorer"):
        total_elapsed = time.perf_counter() - self.start_global
        usage_self = resource.getrusage(resource.RUSAGE_SELF)
        usage_child = resource.getrusage(resource.RUSAGE_CHILDREN)
        total_cpu_time = (
            usage_self.ru_utime + usage_self.ru_stime + 
            usage_child.ru_utime + usage_child.ru_stime
        )

        if self.mini_log:
            msg = f"[Done] {module} completed. (Wall: {total_elapsed:.1f}s, CPU: {total_cpu_time:.1f}s)"
            self.logger.info(msg, extra={'bare': True})
        else:
            msg = [
                f"**** {module} completed! ****",
                f"CPU time  : {total_cpu_time:.1f}s",
                f"Wall time : {total_elapsed:.1f}s"
            ]
            self.logger.info("\n".join(msg), extra={'bare': True})

    def finish(self, msg=f"Done!\n"):
        self.logger.info(msg)

    def warning(self, msg, **kwargs):
        self.logger.warning(msg, **kwargs)

    def error(self, msg, **kwargs):
        self.logger.error(msg, **kwargs)

    def run_cmd(self, cmd, shell=False):
        self._validate_command(cmd)

        if isinstance(cmd, list):
            cmd_str = shlex.join(cmd)
        else:
            cmd_str = cmd

        executable_cmd = cmd_str if shell else cmd
        
        if not self.mini_log:
            self.logger.info(f"Running: {cmd_str}")
        else:
            fh = self.logger.handlers[0]
            fh.emit(logging.LogRecord(self.logger.name, logging.INFO, "", 0, f"CMD: {cmd_str}", None, None))

        try:
            if self.mini_log:
                subprocess.run(executable_cmd, shell=shell, check=True, capture_output=True, text=True)
            else:
                prefix = "  > "
                with subprocess.Popen(
                    executable_cmd,
                    shell=shell, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT, 
                    text=True, 
                    bufsize=1,
                    encoding='utf-8',
                    errors='replace'
                ) as p:
                    for line in p.stdout:
                        msg = line.rstrip()
                        self.logger.info(f"{prefix}{msg}", extra={'bare': True})

                if p.wait() != 0:
                    raise subprocess.CalledProcessError(p.returncode, cmd_str)

        except subprocess.CalledProcessError as e:
            self.logger.error(f"\nError running command: {cmd_str}")
            if hasattr(e, 'stderr') and e.stderr: 
                self.logger.error(f"Details:\n{e.stderr}")
            sys.exit(1)

class FallbackLogger:
    def info(self, msg):
        logging.info(msg)

    def warning(self, msg):
        logging.warning(msg)

    def error(self, msg):
        logging.error(msg)

    def run_cmd(self, cmd, shell=False):
        if isinstance(cmd, str):
            executable = cmd.split()[0]
        else:
            executable = cmd[0]
        if not shutil.which(executable):
            raise ValueError(f"Required command '{executable}' not found in PATH.")

        logging.info(f"Running: {cmd}")
        subprocess.run(cmd, shell=shell, check=True)

    def start_step(self, name):
        logging.info(f"{name}...")

    def end_step(self):
        logging.info(f"Done.")

    def success(self):
        logging.info("Pipeline finished successfully.")

    def finish(self, msg="Done!"):
        logging.info(msg)

    def __getattr__(self, name):
        return lambda *args, **kwargs: None

_fallback_logger = FallbackLogger()

def get_logger(instance=None):
    return instance if instance else _fallback_logger

def init_logger(outdir: Path, prefix: str, debug: bool, mini_log: bool, total_steps: int, version: str, start_step: int = 1):
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    return WorkflowLogger(outdir, prefix, debug, mini_log, total_steps, start_step, version)
