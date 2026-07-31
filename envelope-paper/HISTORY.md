
**Auto-sync gated (decision: gate on the test suite).** The cron fired again at
00:00 and committed this whole round as "Auto-sync: 2026-07-31 00:00", pushed to
GitHub before review — the third demonstration in one session. Already public, so
no history rewrite; 16e1dc7 is this round's work under an unhelpful message.
`git-auto-sync.sh` now runs the 103-test suite in `envelope-paper/` before staging
anything and skips the tick if it fails (verified: breaking a manuscript number
makes it refuse to commit). Two incidental fixes: it no longer logs "Synced" when
`git add -A` stages nothing — which it had been doing every 30 minutes because the
`proteostasis-paper` submodule always reads dirty — and it refuses to commit at all
if the suite cannot be run.
