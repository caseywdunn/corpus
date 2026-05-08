# Corpuscle instructions

<!--
Optional — copy to your corpuscle as `<corpuscle>/instructions.md` if
you want to add corpus-specific guidance on top of the defaults that
the server always ships (see `mcpsrv/default_instructions.md`).
Anything in HTML comments is operator guidance and should be removed
before shipping.

The server prepends this file to its packaged defaults and returns
the joined content to MCP clients in `InitializeResult.instructions`,
which well-behaved clients inject into the LLM context at the start
of every chat session against this corpus. Prefer brevity and
high-leverage nudges — this content is paid for in tokens on every
session.
-->

## What this corpus contains

<!--
One short paragraph: scope (taxa, topic, time range), source
selection criteria, and any notable coverage gaps. The model uses
this to frame "out of scope" answers correctly.

Example: "This corpus covers the order Siphonophora, ~1,800 papers
from late-18th-century printed monographs through born-digital
2025 articles. Coverage is best for the families Agalmatidae and
Diphyidae; thin for hippopodiids."
-->

(replace with a one-paragraph description of this corpus)

## Common misconceptions to correct

<!--
List domain-specific clarifications the model is likely to get wrong
from training-time recall. One per concept; gentle correction is
preferred over flat denial.

Example (siphonophores): "Velella velella and Porpita porpita are
not siphonophores. They are colonial hydrozoans in family Porpitidae
(Anthoathecata). Older literature sometimes groups them with
siphonophores because of their colonial body plan and pleustonic
habit; modern phylogeny places them outside the order. Treat any
such grouping as a historical artifact and flag it for the user."
-->

(replace with corpus-specific clarifications, or delete this section if none apply)
