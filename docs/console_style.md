# Console and logs

What rpsHunter writes to the terminal and to its logs follows one convention. Its executable form
is `workflow/scripts/log.py` (what jobs write), `workflow/scripts/external.py` (running tools) and
`workflow/scripts/console.py` (what the launcher shows).

## Rules

1. **Jobs write plain lines; only the launcher draws.** A script never colours anything.
2. **One line, one shape.**

   ```
   14:02:40 OK select_hits Human_PRDM9/Homo_sapiens | 33,721 of 33,863 hits selected
   time     level step     scope                     message
   ```

   The scope is `<query>/<species>`, or `all` for a job over the whole run. With many jobs writing
   at once, step and scope are what make interleaved lines readable.
3. **A level means one thing.**
   - `ERROR`: the job cannot produce a correct output. It stops with a non-zero status.
   - `WARN`: the output exists but a person should look. A clean run of the model genomes shows no
     warning, so a warning always deserves a look.
   - `OK`: the one headline a job prints when it finishes, with its key number.
   - `INFO`: detail, kept in the job's log.
4. **Every line is kept.** The screen shows what the verbosity asks for; the logs hold everything.

## Where things go

| What | Where |
|---|---|
| A job's lines and its tool's error output | `logs/jobs/<rule>/<query>/<species>.log`, one file per job |
| Everything a run printed, ours and Snakemake's | `logs/runs/<date>_<time>.log` |
| The warnings of a run, when there are any | `logs/runs/<date>_<time>.warnings.txt` |
| Snakemake's benchmark of each job | `logs/benchmarks/<rule>/...tsv` |

## The screen

```
rpsHunter 0.9.0 | config /path/to/config.yaml | 16 cores | --search --domains --loci
✓ tblastn  Human_PRDM9/Homo_sapiens  33,863 hits at E <= 10
✓ select_hits  Human_PRDM9/Homo_sapiens  33,721 of 33,863 hits selected
...
Done in 11 min 32 s: 78 of 78 jobs, no warnings.
Run log: logs/runs/2026-11-27_004918.log
```

| Verbosity | Shown |
|---|---|
| `quiet` | errors, warnings, the progress line, the summary |
| `normal` | also one line per finished job |
| `verbose` | every line, Snakemake's included |

Set it with `display.verbosity` in the config or `--verbosity` for one run. On a terminal a
progress line (`[12/78] 15%`) stays below the messages.

When a job fails, its `ERROR` line is shown once, and the summary names each failed rule with the
path of its log, where the tool's own messages and the traceback are.

## Colour

Colour carries the level and nothing else, and never alone: each level also has a symbol, or a
word where the terminal cannot show symbols. Colours are the terminal's own, so they follow its
theme. There is no colour when the output is not a terminal; `NO_COLOR=1` and `FORCE_COLOR=1` turn
it off and on. Files never contain colour codes.

## Words

- Say what happened in the reader's terms, with the number that matters.
- Counts carry thousands separators: `33,721`.
- A known failure says what to do: `PipelineError("what went wrong", hint="what to do")` shows as
  one line ending in `Fix: ...`.

## Adding a script

1. `log.add_scope_argument(parser)` before parsing, `log.job_logging("<step>", args.scope)` after.
2. End `main` with one `log.ok(...)` headline; run it with `sys.exit(log.run_main(main))`.
3. In its rule: `log: job_log(...)`, pass `--scope`, and end the command with
   `2>&1 | tee {log:q}` so the lines reach both the job's log and the launcher.
4. Run tools through `external.run_tool`: a failure then stops the job with one line.
