# Study questions and decisions subsystem design

> **Source implementation and longer-term design.** The deterministic first
> tranche described below is implemented in the source build. Subsequent design
> sections specify the broader target system; their additional commands,
> configuration fields, and schemas are proposals, not current CLI contracts.
> No published-release or remote LLM execution claim is made here.

JLinAlg should resolve what it can from evidence, ask the scientist only for
decisions or knowledge it cannot safely obtain, and keep independent work moving.
Built-in rules and an optional LLM can propose questions; the same deterministic
validation, permission, amendment, and provenance rules govern both. An answer
is not automatically a fact, an approval, a scientific conclusion, or permission
to execute code.

This design extends the [auditable study kernel](auditable-omics-pipeline.md)
and its [tutorial](study-pipeline-tutorial.md). It applies to any supported assay,
tissue, phenotype, and starting stage, not to a particular demonstration study.

## Implemented first tranche

The source build documents scientist questions before making them available for
review. It does not claim a question was actually read. Questions have a stable
semantic ID, immutable revision, exact prompt, reason, typed response contract,
affected hypotheses/operations/datasets, input hashes, and originating rule.
Runtime questions also point to hashed evidence manifests. Unknowns and optional
work that was declined remain distinguishable from negative scientific evidence.

The first built-in rules cover missing phenotype files/outcomes for association,
shared-parent cohort relationships, declared dependence that the built-in pooling
adapters cannot handle, missing or ambiguous feature annotation for genetic
follow-up, GWAS/QTL input specifications, and optional diagnostic review. They
use explicit fields and typed `StudyNeed` records, not guesses from error text.
They are not a general intake solver: malformed/incomplete core contracts still
fail validation, and some older blockers still appear only as operation requests.

`study questions` reads the local inbox; `--out answers.yaml` exports blank
responses without preaccepting anything or overwriting an existing file.
`study answer answers.yaml --actor-id PERSON --check` validates the complete
batch without writing. Replace `--check` with `--record` to preserve the original
submission, normalized answers, rationale, citations/source notes, and optional
small local metadata snapshots. Local-file answers must refer to existing files
inside the project; evidence snapshots are UTF-8 metadata of at most 1 MiB.

An answer does not execute an analysis, grant disclosure permission, or change
the resolved study. Supported factual changes produce one proposed
`amendments/proposed/amendment-ID.yaml` for the batch; supply a user-chosen
`amendment_id`. Review with `study amend ... --check`, then apply with
`--approved-by PERSON`. Plan approval and execution follow separately (or through
the existing `study auto` command). Narrative review answers require a separately
prepared amendment or supported inputs; no free text is interpreted as code.

The currently supported proposal targets are dataset `files.phenotype`,
`files.annotation`, `genome_build`, `coordinate_system`, and `units`; hypothesis
`outcome`; and explicit pairwise cohort relationships. Not all targets have a
built-in question trigger yet. Broader hypothesis and analytical-plan changes
continue to use the existing amendment format.

Questions and answers are rendered in `QUESTIONS.md`, `SESSION_HISTORY.md`,
`state/questions.yaml`, project-status counts, and the searchable HTML/JSON index.
The event ledger anchors their immutable artifacts and actor/timestamp records.
The local index includes sensitive answer text; it is not a disclosure bundle.
Raw answers and evidence snapshots are not added to `index/LLM_CONTEXT.md` or
to the existing allowlisted assistance previews. No LLM is called by the question
coordinator or answer importer.

The optional advisor operation can now propose up to three labelled scientist
questions, with provider/model and interaction provenance, for this same local
inbox. These are review questions without amendment targets or execution authority.
Related phenotype/source suggestions and hashed code text also remain proposals;
see the [bioinformatics workflow](bioinformatics-evidence.md#questions-and-optional-ai).

The executor gates affected operations, not the entire project. Accepting general
plan assumptions does not answer a shared-parent cohort question or make related
cohorts statistically independent. `unknown` and `defer` retain the prerequisite;
`decline` is available only for optional review. None is treated as successful
completion. Repeating a run retains unanswered/answered questions when their
premises match. Ordinary answer exports contain only new open questions;
`--include-answered` also exports unresolved answered/deferred/declined questions
for explicit correction, naming the previous answer in `supersedes`.

Changed inputs or study contracts invalidate old responses. A correction needs
a new answer ID and explicit `supersedes`; the previous answer remains visible.
An identical repeated submission is idempotent, while an ID reused for different
content or attribution is rejected. Whole-study cache invalidation after an
amendment remains conservative; selective scientific dependency invalidation is
not implemented.

Project mutations share a lock. Answer batches are validated before publishing
immutable files and one ledger commit; uncommitted files grant no authority.
An identical retry can finish interrupted publication. `study verify --repair-views`
rebuilds derived notebooks/status/indexes only after authoritative ledger records,
question/answer artifacts, and analysis evidence verify. It cannot repair missing
or altered authoritative evidence. A general crash-recovery protocol for all
older plan/amendment/run operations remains future work.

These are local attribution and integrity safeguards, not authentication,
digital signatures, regulatory certification, or a sandbox for external code.
An administrator could rewrite an unsigned history. Authenticated remote
approval, interactive prompting, LLM-generated question proposals, intake of
incomplete contracts, and generated-code execution remain later milestones.

See the [question workflow tutorial](study-pipeline-tutorial.md#documented-scientist-questions)
for the actual answer format and commands. The sections below retain the broader
target design; where they differ, this implemented contract takes precedence.

## Design commitments

- Ask a question only when its answer could change an admissible action,
  interpretation, or explicitly requested research objective.
- Separate factual clarification, scientific investigation, and authorization.
- Work without an LLM. LLM assistance may broaden proposals but cannot weaken
  scientific prerequisites or acquire approval authority.
- Keep study descriptions small. Generate question records and implementation
  detail rather than asking the researcher to author them.
- Bind decisions to the evidence and scope actually reviewed. Never silently
  reuse an answer after its relevant premises change.
- Preserve unknowns, negative evidence, declined work, and changed hypotheses.
- Never interpret silence, an expired prompt, a recommendation, or `mode: auto`
  as consent to a new disclosure, hypothesis, expense, or code execution.

## What the current implementation provides

The following integration points were checked in the local source on
2026-10-01. This is an implementation inventory, not a claim about a published
release.

| Existing component | Available foundation | Remaining extension |
| --- | --- | --- |
| `StudyValidator` and `StudyResolver` | Validated contracts, dictionary mappings, explicit defaults and relationship assumptions | Typed issues and an intake path for incomplete but parseable contracts |
| `StudyPlanner` | Immutable plans, assumption approval, deterministic prerequisite questions | Broader scientific question rules and explicit decision dependencies |
| `StudyOperationExecutor.Result` | Complete, not applicable, waiting, failed, and typed needs | Migrate remaining untyped blockers |
| `StudyObservations` | Immutable observations and proposed amendments | Observation-to-question links and competing explanation records |
| `StudyAmendments` | Base-hash validation, impact preview, answer-linked proposals and immutable applied amendments | Broader typed decision proposals |
| `StudyLedger` and `StudyProject` | Hash-chained events, questions/answers, actor attribution, derived notebooks, shared mutation lock | General transaction recovery and authenticated authority |
| `StudyExecutionEngine` | Scoped prerequisite gates, verified resume, immutable attempts and input/implementation hashes | Selective decision-aware reuse after amendments |
| `StudyAssistance` | Opt-in bounded provider requests with request/response/decision files | A separate, explicitly enabled question-proposal capability |
| `StudyEvidenceIndex` and `StudyCli` | Local HTML/JSON inbox, blank answer export, checked batch import | Interactive presentation and remote approval channels |

Current limitations matter: legacy `pending_decisions` is separate from the new durable question inbox;
`--approved-by` records a supplied name rather than authenticating a person; and
the current execution cache includes the whole resolved-study hash. Selective
reuse after an amendment therefore requires additional dependency work. Current
LLM adapters do not execute generated code.

## Route each question to the right respondent

A research question is not necessarily a prompt to a human. The coordinator
routes a typed need through the following destinations.

| Need | First destination | Human question if still needed |
| --- | --- | --- |
| Missing assembly, effect units, phenotype meaning, or cohort relationship | Approved contract, dictionary, source metadata and previous applicable answers | Ask for the missing fact and its source; permit an explicit unknown |
| Ambiguous or conflicting metadata | Local validation and a conflict report | Ask which source is correct; do not silently choose a default |
| Possible batch effect, cell-composition effect, or unusual pattern | Eligible diagnostic or sensitivity analysis | Ask for unavailable inputs or permission for a new analysis, not whether the scientist thinks the result is an artefact |
| Which external tissue or phenotype could inform a hypothesis | Registered relevance rules, verified literature, optional LLM proposal | Ask about unresolved scope or competing scientific objectives |
| A new hypothesis or changed primary analysis | Amendment proposal and impact preview | Ask the scientist to approve, revise, decline, or defer the change |
| Remote transmission, resource cost, or an untrusted executable extension | Effective permission policy and execution capability checks | Ask only when a permitted approval path exists; a denial is not a reason to solicit repeated exceptions |
| Network failure, missing credential, unavailable service | Bounded retry policy and actionable status | Request a local resource or credential setup when useful; never request a secret in an answer file |
| Explanation of an existing result | Verified local evidence and optional bounded explanation | Present an explanation, not a blocking question |

Approved inputs and accepted amendments take precedence over suggestions. A
data dictionary can resolve a missing mapping, but a contradiction with the
approved contract becomes a conflict. A guessed assembly or an LLM confidence
score cannot supply a missing scientific prerequisite.

Existing default independence for undeclared cohort relationships is retained
as an explicit assumption where applicable. Evidence of overlap, a shared
parent cohort, or a declared relationship triggers review of affected joint
analyses. Choosing to proceed under an assumption never relabels it as verified
independence, and cannot bypass a method's unmet requirements.

## How questions are generated and selected

The proposed processing sequence is:

~~~text
contract or runtime issue, observation, or scientist request
  -> local fact lookup and applicable prior decisions
  -> built-in rule candidates plus optional LLM candidates
  -> schema, evidence, relevance, feasibility and policy checks
  -> deduplication, grouping and priority
  -> tool investigation OR scientist question OR recorded disposition
  -> answer validation and impact preview
  -> required amendment and approval
  -> replan and resume eligible work
~~~

### Built in rule registry

Each versioned rule declares its trigger code, supported inputs, relevant
scientific scope, fact dependencies, response schema, permitted action intents,
authority requirement, blocking scope, and test fixtures. Rules are data-driven
where possible; modality-specific plugins contribute rules through the same
interface. Example initial rule families are:

- Ambiguous phenotype or effect scale before harmonization or pooling.
- Missing assembly or feature annotation before coordinate-dependent follow-up.
- Known or suspected cohort dependence before independent-cohort inference.
- Missing tested-feature universe before a background-dependent enrichment.
- Missing QTL or regional GWAS inputs before eligible genetic follow-up.
- Diagnostic observations that warrant an available sensitivity analysis.
- An unsupported requested operation or missing external provenance.
- A hypothesis whose prerequisite is no longer supported by available evidence.

Rules must not make all follow-ups prerequisites of the main analysis. For
example, a missing protein mapping should not block valid RNA-level analyses.

### Optional LLM proposal provider

The advisor receives a policy-approved context bundle and the allowed proposal
schema. It can suggest a question, competing explanations, relevant sources,
candidate tissues, and a possible analysis. It cannot declare an answer correct,
set an approval requirement, disable a gate, or supply an executable command as
an action. The engine assigns those properties independently.

Every proposal must identify a hypothesis or a project-wide dependency, explain
how an answer could change the next step, distinguish observations from
speculation, and cite the evidence it relies on. External citations are resolved
and checked; a missing or irrelevant citation remains unverified. A creative
proposal without prior literature may be retained as an explicitly speculative
new hypothesis, not presented as an established mechanism.

Retrieved pages, papers, dictionary prose, and provider responses are untrusted
content. Instructions embedded in them cannot alter the study, authority,
disclosure policy, or command allowlist. Escaped rendering and bounded parsing
are required for every presentation surface.

The engine preserves accepted, rejected, duplicate, and unsupported proposals
with reasons. A question never triggers repeated LLM calls on unchanged context.
An unavailable LLM does not prevent deterministic questions or local execution.

The approved advisor policy must also bound calls, tokens, elapsed time, and
cost where applicable. Record the effective limits and consumption; exhausting
a limit stops optional advice, not eligible deterministic work. If a required
cost or execution bound cannot be enforced, retain the proposal locally for
review rather than starting an unbounded call. The three-question display limit
below is not a substitute for these execution budgets.

### Reduce interruption without concealing uncertainty

Before showing a question, the coordinator asks whether the answer is already
available, whether an authorized local check can obtain it, whether another
question must be answered first, and whether the answer would actually affect
the current plan. It collapses repeated downstream needs into one root question.

Use explainable priority bands: safety and contradictory premises; missing facts
blocking requested analyses; approvals for requested work; then optional
exploration. Within a band, prefer questions that unblock more independent work,
then stable identifiers. Do not invent numerical information-gain probabilities.

Present up to three questions per interactive round by default. This is a
display limit, not a limit on retained blockers. Include the full queue in the
index. Related questions may share a display group but retain separate answers
and permissions. Never bundle remote-data consent into an ordinary scientific
choice. Deferral suppresses repeat prompts until its declared revisit condition;
it does not grant permission or erase the blocker.

## Question and answer contracts

Use YAML/JSON records with a version independent of `study.yaml`. The user
normally sees a short question card and edits an exported answer template;
JLinAlg generates identifiers, hashes, scope, and provenance.

### Question record

| Field group | Required meaning |
| --- | --- |
| Identity | Project ID, stable question ID, revision, semantic key, record hash and schema version |
| Origin | Rule ID/version or LLM interaction reference, proposer actor/model, creation event and UTC timestamp |
| Scientific scope | Hypothesis IDs, dataset/feature IDs where applicable, and primary, sensitivity, exploratory, or alternative-explanation role |
| Premises | Resolved-study/plan snapshot, relevant dependency hashes, evidence references with row or field selectors, and verified source citations |
| Presentation | Plain-language question, why it matters, known facts, uncertainty, and consequences of each option |
| Response | Typed answer schema, stable option IDs, explicit unknown/defer handling, and optional explanation or resource reference |
| Decision boundary | Engine-assigned authority requirement, affected operations, action intents, and what can continue without an answer |
| Recommendation | Optional recommendation with rationale and origin; never a preaccepted answer |

A semantic key identifies the issue, such as a dataset's assembly, rather than
its wording or provider. A revision pins the exact premises and consequences.
Presentation-only paraphrases cannot change the machine contract. Deduplication
uses normalized scientific scope and decision target, not an LLM similarity
judgment alone.

The following abbreviated generated example uses symbolic hashes and operation
IDs. It is illustrative, not a valid current input file.

~~~yaml
schema_version: 1
id: Q0001
revision: 1
kind: fact_clarification
semantic_key: dataset:d1:genome_assembly
origin:
  type: rule
  id: annotation.assembly_required
  version: 1
scope:
  hypotheses: [H1]
  datasets: [d1]
  role: prerequisite
context:
  resolved_study_hash: <hash>
  dependency_hash: <relevant-premises-hash>
evidence:
  - artifact: validation/d1-annotation.json
    sha256: <artifact-hash>
    selector: /genome_assembly
prompt: Which genome assembly was used for the coordinates in dataset d1?
why: Regional QTL and GWAS lookups require an explicit compatible assembly.
response:
  type: single_choice
  options:
    - id: GRCh37
      label: GRCh37
    - id: GRCh38
      label: GRCh38
    - id: other
      label: Another assembly with supporting metadata
    - id: unknown
      label: I do not know
  allow_defer: true
  source_required_for_known_value: true
recommendation: null
gate:
  operations: [regional-qtl-lookup, regional-gwas-lookup]
  otherwise: waiting_for_input
  independent_work: [association, replication, meta_analysis]
resolution:
  intent: propose_study_amendment
  authority: scientist
~~~

Single choice, multiple choice, quantities with units, structured resource
references, and constrained text are supported response shapes in the design.
Dates, allowed values, ranges, and cross-field conditions are checked explicitly.
Free text is not executable and cannot silently become a structured instruction.
If an LLM interprets free text or a meeting note, show the proposed structured
answer and its consequences for confirmation. Preserve the original wording.

### Answer and decision records

An answer binds to a project, question revision and question hash. It records
the value, explanation, source references, submitting actor, claimed respondent,
server-recorded UTC time, and optional meeting ID. Actor authority comes from
the trusted invocation or authenticated session, not a field supplied by an LLM.

A minimal exported answer template would look like this:

~~~yaml
schema_version: 1
project_id: example-study
responses:
  - question_id: Q0001
    question_revision: 1
    question_hash: <hash-from-export>
    disposition: answer
    value: GRCh38
    rationale: Confirmed in the original processing manifest.
    evidence:
      - file: metadata/original-processing.yaml
        selector: /reference/assembly
~~~

File references are resolved and hashed locally under an approved read scope;
they are not fetched remotely or passed to an LLM merely because they occur in
an answer. Rejected answers receive validation feedback without mutating the
scientific contract. An accepted scientist assertion remains labeled as such if
it cannot be independently verified.

The derived decision record separately contains the normalized answer, its
validation result, effective policy, proposed action/diff, approval references,
and final disposition. Permission questions require a separate explicit approval
of the exact action and scope. Supplying a fact does not grant remote access.

## Lifecycle and scientific changes

Question state and operation state are separate. An answer can be valid while
an affected operation still lacks another prerequisite.

| Question state | Meaning and allowed next step |
| --- | --- |
| `open` | The need is current; collect a response or resolve it from verified evidence |
| `answered` | A response was recorded; validate its type, premises, and authority |
| `needs_clarification` | The response is incomplete, conflicting, or cannot be safely interpreted |
| `awaiting_approval` | The response implies a concrete amendment or controlled action awaiting approval |
| `resolved` | A verified fact, applicable policy decision, or approved change addresses the question |
| `deferred` | Still unresolved; reconsider only on its recorded condition or explicit scientist request |
| `declined` | A proposed optional action was explicitly rejected; no action is authorized |
| `superseded` | Premises or decision scope changed; retain history and link any successor |
| `withdrawn` | The question is no longer relevant, with a recorded reason |

An explicit `unknown` answer records useful information but does not satisfy a
required assembly, unit, or relationship fact. Record it, retain the gate, and
offer a feasible alternative or defer further prompting. A recommendation and
its selection probability are never substitutes for an answer.

Changing the study contract produces a proposed `amendment-XXX.yaml` with a
user-selected ID. The proposal names the originating question, answer, evidence,
and rationale, and previews affected hypotheses, analyses, assumptions, and
results. The original hypothesis and earlier findings remain immutable. Adding
a sensitivity analysis, replacing an unsupported premise, and creating a new
exploratory hypothesis have distinct dispositions.

The existing amendment validation and plan approval remain authoritative. No
answer directly edits `study.yaml`, `resolved-study.yaml`, or a running plan.
An explicit authorization to apply a reviewed amendment is not permission to
run an arbitrary extension, purchase a service, or disclose data. These have
separate policy gates.

An accepted fact already present in the approved contract can close a redundant
question without an amendment. Annotations and explanations that do not alter
analysis semantics can remain ledger events. Any change to dataset interpretation,
scientific scope, execution policy, or analysis configuration must go through
the appropriate amendment and approval path.

### Corrections and concurrent answers

Accepted answers are immutable. Corrections create superseding answer and
decision records with a reason and impact preview; prior results are retained
and marked as based on the superseded decision. Revoke pending permissions
before additional affected work starts; already transmitted data cannot be
recalled by changing an answer.

Import batches atomically after cross-answer validation. Duplicate submissions
with the same idempotency key and content return their existing receipt. Reuse
of that key with different content is rejected. Two incompatible answers do not
use last-writer-wins: return a conflict and require explicit reconciliation.

Bind approval to the exact question revision, action proposal, policy, and
relevant evidence hashes. A stale answer cannot approve a revised question.
Unrelated project events need not invalidate an answer; the engine checks its
actual dependency set. Amendments still use the current full base-study hash.

## Execution gates and resume

At validation, planning, and operation boundaries, the coordinator reconciles
needs and decisions. An operation proceeds only when its scientific inputs and
applicable permissions are satisfied. A question may block one operation, a
dependent branch, or, for a truly global safety/integrity failure, the whole run.
Independent work and evidence indexing continue whenever safe.

Batch mode never waits on standard input. It records the inbox, executes all
eligible work, writes status/index/notebook views, and exits. Interactive mode
asks a short batch at a safe checkpoint, with a cancel/defer option. Responses
cannot mutate an in-flight plan; a changed plan starts at the next execution
boundary. The initial implementation can collect answers after the run rather
than adding a background service.

On a repeated command:

1. Verify the ledger and all referenced decision artifacts.
2. Reconcile questions against current inputs, accepted answers, and policies.
3. Preserve unanswered and deferred questions without duplicate prompts or
   repeated unchanged LLM proposals.
4. Apply no unapproved change; load only an approved current plan.
5. Reuse verified eligible work and retry affected unfinished operations.
6. Report remaining questions and precisely which requested objectives are unmet.

For the first implementation, retain the current conservative cache behavior:
an amendment may invalidate more work than strictly necessary. A later selective
cache must include all relevant parameters, artifact content, decisions, policy
effects, and implementation identity. Changes to gate authority must be checked
even when a numerical artifact is reusable. Reused artifacts keep their original
provenance and acquire an explicit new-plan reuse record.

`--rerun` means recompute; it does not mean forget answers, change hypotheses,
erase history, or grant permission. An optional action declined by the scientist
is `not_pursued`, not scientifically `not_applicable` and not a completed analysis.
That disposition needs a new scheduler/report concept; it must not be faked
using the current executor's `NOT_APPLICABLE` state. No required scientific
prerequisite can be bypassed by declining its question.

Transient network errors, policy denial, credentials missing, successful empty
search, missing source coverage, incompatible data, and lack of an implementation
are distinct typed reasons. Repeating an offline command must not repeatedly ask
the researcher to override an intentional denial.

## Human interfaces and optional configuration

### Question card

The HTML index and interactive CLI present the same content:

- The question and its hypothesis or project-wide scope.
- Why it is being asked now, with links to supporting and conflicting evidence.
- Options and their scientific, computational, privacy, and cost consequences.
- A justified recommendation when appropriate, clearly attributed.
- What stops without an answer and what can continue.
- How to answer, supply a source, say unknown, decline an optional proposal, or defer.

The index shows open, deferred, awaiting-approval, resolved, and historical views.
Search filters include hypothesis, dataset, question type, source, actor/model,
and time. The notebook links every question, answer, approval and resulting run.
No browser-side script can silently apply an answer to local files. Initially,
the HTML view is read-only and exports an answer template for CLI submission.

### Proposed command sequence

These are proposed commands, not runnable instructions for the current JAR:

~~~text
study questions --project PROJECT
study questions --project PROJECT --out answers.yaml
study answer answers.yaml --project PROJECT --check
study answer answers.yaml --project PROJECT --record --actor-id scientist
~~~

`questions` reads the derived inbox and optionally exports a blank template; it
does not contact providers or generate new questions. `answer --check` previews
validation and consequences without recording an answer. `answer --record`
records validated responses and creates immutable proposals where needed; it
does not execute analyses or grant permission. Repeated exports preserve stable
IDs/hashes, and refuse to overwrite an existing answer file without an explicit
file-replacement option.

If an answer changes the contract, the subsequent workflow uses the existing
`study amend --check`, `study amend --approved-by`, `study plan`, and
`study approve` gates before the next `study run` or `study auto`. The new question
layer must not let `--accept-assumptions` approve arbitrary question types.
Granular remote-request and extension approval interfaces belong to later
milestones, not to the first question CLI.

Proposed configuration is intentionally small; it is not accepted by the current
version-1 validator:

~~~yaml
interaction:
  mode: batch
  max_questions_per_round: 3
assistance:
  questions: false
~~~

Omission would select these defaults in the future schema. Existing assistance
fields still control provider configuration. Enabling `assistance.questions`
also requires assistance itself to be enabled and the effective disclosure
policy to permit the request. An upgrade must not silently expand an existing
LLM permission to this new purpose. Interactive presentation is opt-in, for
example through a proposed `study auto --interactive` override, and requires a
terminal; unattended environments fall back to batch.

For question-enabled execution, the proposed exit contract is 0 when requested
work is complete or only explicitly dispositioned optional work remains, 3 when
requested work awaits input/approval, and 2 for an execution or integrity failure.
Machine output must retain complete counts and reasons regardless of exit code.
Today `study run` can return 3 for waiting work, while `study auto` can return 0
despite waiting. Unifying those semantics is a documented, versioned migration,
not a silent change in this design-only work. Optional exploration proposals
alone should not turn a completed requested analysis into a failed run.

## Provenance and security boundaries

Use the existing event ledger as the authority. New immutable artifact folders
hold questions, answers, decision proposals and reviewed decisions. For example:

~~~text
questions/Q0001/revision-0001.yaml
answers/A0001.yaml
decisions/D0001/proposal.yaml
decisions/D0001/resolution.yaml
amendments/proposed/amendment-XXX.yaml
ledger/events.jsonl
state/questions.yaml
state/project-status.yaml
SESSION_HISTORY.md
index/index.json
index/index.html
index/LLM_CONTEXT.md
~~~

`state/questions.yaml` is a rebuildable view, not a second source of truth. Event
types cover question creation/supersession, answer recording/validation, decision
proposal/approval/decline, amendment linkage, and resulting execution. Resolution
artifacts are created once, not edited as status changes. Full model interaction
records remain separate and linked by hash, with requested/returned model IDs,
timestamps, prompts, replies, disclosure manifests, and accepted/rejected
proposals. Human-readable decision rationale is required; private model chain
of thought is neither needed nor represented as available.

The entity, activity, and agent relationships can map to
[W3C PROV-O](https://www.w3.org/TR/prov-o/). This is a proposed interchange mapping,
not a claim of implemented PROV export or regulatory certification.

All mutations need one project-wide transaction protocol, including answer
imports, amendments and plan changes, not only the existing run/ledger locks.
Stage immutable files; validate their hashes and expected current revisions
under the shared lock; publish the files; append a commit event referencing
them; then regenerate disposable views. Uncommitted files left by a crash grant
no authority. Recovery identifies them and either resumes the same idempotent
transaction or quarantines them. A missing file referenced by a committed event
is an integrity failure. A long-running execution must not hold a lock while
waiting for a person's answer; the first tranche records answers between runs.

Local actor names provide attribution, not proof of identity. The first release
must state that limitation. An external LLM, MCP client, or unattended service
must not gain human authority by sending `actor_type: human` or a name string;
remote approval requires an authenticated, capability-scoped channel in a later
milestone. Local filesystem administrators remain part of the trust boundary.
Hashes detect changes relative to a trusted ledger/head; they do not prevent an
administrator from rewriting an entire unsigned history.

Questions and answers can themselves contain sensitive research information.
Keep the full notebook and raw answers local by default; do not automatically
copy them into an externally eligible context bundle. A local LLM-oriented index
is not itself permission to transmit. Apply denial precedence and purpose-bound
egress policy to all providers and search services. Never log API keys; require
approved file roots, bounded document parsing, escaped HTML, and secret checks
on outgoing context. Exported evidence packages need a separate disclosure and
licensing review; a locally complete archive is not automatically shareable.

## Scientific safeguards and examples

### Whole blood expression and sleep

For an illustrative H1 linking whole-blood expression to sleep, the engine
first resolves what sleep means in the study: duration, insomnia, chronotype,
or another measured phenotype. It uses the dictionary when authoritative and
asks only if the distinction remains unclear. Related phenotypes are not
silently substituted for the stated outcome.

An investigation proposal could include blood, immune-cell, and brain QTL
contexts with distinct purposes. Brain follow-up is mechanistic exploration,
not replication of the original blood association. A source record must include
its original study citation, resource/reprocessing citation, accession/release,
relevant tissue and phenotype, sample and ancestry information, data limitations,
and a hypothesis-linked selection rationale. Selection precedes viewing the
desired association result where feasible. The eQTL Catalogue explicitly asks
users to cite both the resource and the original contributing studies
([citation guidance](https://www.ebi.ac.uk/eqtl/Studies/)).

If cross-tissue exploration is already authorized within the approved scope,
eligible registered steps can proceed without a redundant question. Otherwise,
ask whether to add a precisely described exploratory branch, showing candidate
contexts, evidence, tests, and expected costs. Lack of a significant QTL is not
converted into a claim that a biological mechanism is absent; report coverage,
power limitations where assessable, and the analysis actually performed.

### An unusual association pattern

An observation records the pattern and exact artifact. The system proposes
competing biological and technical explanations. It runs already approved
diagnostics where possible. If a batch explanation needs metadata that are not
available, it asks for the relevant manifest or records an unresolved limitation.
It never asks the scientist to certify that a pattern is biological merely to
unblock a preferred narrative. Summary-only inputs remain summary-only.

### A failed prerequisite hypothesis

If H2 assumes evidence from H1 that the current results do not establish, retain
H1's results and the original H2. Present options to retain H2 as unsupported,
retire it with a reason, or propose a revised/new hypothesis. Do not equate a
nonsignificant result with proof that H1 is false. New objectives receive
explicit hypothesis IDs, origin links, and exploratory status as appropriate.

### A proposed new bioinformatics script

An LLM may propose an investigation outside the operation registry. Initially,
record it as unsupported with its input and implementation requirements. A later
extension subsystem can request code, preserve its prompt/code/dependency hashes,
validate it, and obtain separate execution approval for an isolated environment.
Approval of the scientific question does not approve the script. If an adequate
sandbox is unavailable, execution remains blocked; it does not fall back to an
unrestricted shell. Preserved reviewed code supports replay without a new LLM
call. Tests alone do not establish scientific validity.

### Statistical accountability

Record which results were visible when a question, hypothesis, or tissue/trait
search scope was proposed. Keep attempted analyses and negative results, not
only successful follow-ups. Planned and exploratory test families and their
multiplicity strategy must remain distinguishable. Retrospective approval does
not turn a post-hoc hypothesis into a prespecified one, and an LLM cannot select
the testing family after seeing which correction gives significance. Reused
cohorts or the same source mirrored by multiple repositories are not independent
corroboration. The design does not claim generic statistical error control for
arbitrarily adaptive exploration.

## Proposed components and delivery stages

These names are proposed implementation boundaries, not existing public APIs.

| Component | Responsibility |
| --- | --- |
| `StudyNeed` and `StudyQuestion` | Typed missing inputs, questions, immutable revisions, response contracts and evidence references |
| `StudyQuestionRuleRegistry` | Versioned built-in rule descriptors and deterministic trigger evaluation |
| `StudyQuestionAdvisor` | Optional constrained LLM proposal generation with existing disclosure controls |
| `StudyQuestionCoordinator` | Fact lookup, routing, deduplication, dependency-aware priority and question reconciliation |
| `StudyAnswers` and `StudyDecisionPolicy` | Answer validation, conflict handling, authority checks and consequence preview |
| `StudyQuestionStore` | Ledger-anchored artifacts, reducers, recovery and verification |
| `StudyQuestionPresenter` | CLI/HTML/JSON projections sharing the same validated question model |

### First tranche without an LLM

Implement typed needs, rule-based questions, immutable answers, amendment links,
batch CLI, index/notebook integration, and branch-specific gates. Start with
assembly, units, phenotype mapping, cohort dependence, missing analysis inputs,
and unsupported requested operations. Test at least a matrix-first study, a
result-first study, and a variant-based study. Keep current conservative cache
invalidation; do not include generated-code execution.

Intake must distinguish a parseable incomplete contract from corrupt YAML or an
unsafe file. A valid project can use its ledger immediately. Before the first
successful resolution, store an intake snapshot and question records bound to
its raw-file hash in a separate draft area; these are not an approved study or
execution state. Editing that draft supersedes stale intake questions. Import
its history into the project ledger when a validated base contract is resolved.
No analysis, provider call, or remote metadata lookup runs merely to repair an
unapproved draft. The existing rule that changing a resolved base study requires
an amendment remains intact.

### LLM assisted questions

Add the explicitly enabled advisor, citation/source checking, free-text
interpretation with confirmation, scientific relevance proposals, and bounded
question budgets. Record one proposal attempt per relevant context unless the
scientist requests reconsideration or a recorded retry condition changes.
Evaluate against both useful suggestions and adversarial proposals.

### Broader investigation and integrations

Add approved source-acquisition and cross-tissue investigation plugins, stronger
authenticated approval channels, selective cache invalidation, and reviewed
generated-code extensions. CLI, a future GUI, and MCP clients share the same
question/decision API; none bypasses its authority and provenance checks.

## Acceptance tests before implementation is considered complete

1. An authoritative dictionary resolves a mapping without asking; a conflicting
   dictionary produces a question without changing the approved contract.
2. Missing assembly blocks regional follow-up but preserves eligible association,
   replication, and meta-analysis work.
3. Unknown, deferred, declined, unanswered, and approved have distinct effects;
   timeout or EOF never approves an action.
4. Repeating a command creates no duplicate question on unchanged premises,
   does not re-ask a deferred question, and does not repeat unchanged AI calls.
5. Changed relevant inputs supersede the question and invalidate stale answers;
   unrelated events do not discard valid answers.
6. Answer import validates the entire batch, detects cross-answer conflicts,
   rejects stale revisions, and handles idempotent duplicate submissions.
7. A contract-changing answer cannot execute until its amendment and required
   plan approvals are recorded; blanket assumption acceptance cannot bypass it.
8. Concurrent submissions and injected crashes before/after file publication,
   ledger commit, and view generation preserve a recoverable, verifiable state.
9. Generated indexes and the notebook link every decision to its actor, evidence,
   answer, amendment, execution and outputs; tampering is detected.
10. A declined optional investigation is not marked biologically inapplicable;
    unmet requested objectives remain visible with accurate exit semantics.
11. LLM-off and deny-all modes work without provider calls. Missing credentials
    and blocked networks are distinguished from successful empty searches.
12. Prompt injection, forged human actors, secret-bearing answers, unsafe paths,
    malicious HTML, oversized documents, and executable action strings cannot
    cross the authority or disclosure boundary.
13. A brain QTL follow-up is not labeled blood replication; source mirrors and
    overlapping cohorts are not counted as independent support.
14. Unsupported code proposals cannot launch a process; replay of a later
    approved extension uses the recorded code rather than regenerated code.
15. Old project ledgers remain verifiable; no historic human response, citation,
    authentication, or approval is fabricated during migration.

Question quality evaluation should also measure unnecessary interruptions,
missed scientific blockers, usefulness of explanations, and correct attribution.
Fewer questions alone is not a success metric. A smaller reviewed test corpus
across distinct study types is preferable to an unsupported coverage claim.

## Recommended design choices for review

Adopt batch-first interaction, a deterministic question kernel, typed answers,
explicit unknowns, and immutable decision provenance. Preserve current defaults
unless contradicted by evidence, and make LLM question generation a separate
opt-in capability. Defer authenticated remote approvals and generated-code
execution until their safety and scientific validation requirements are met.

The decisions to confirm before implementation are whether these interaction
defaults are suitable, which built-in question families belong in the first
tranche, and which existing commands should opt into the new waiting exit
semantics. The proposed first tranche above is the recommended starting scope.

<!-- SCIENTIFIC-CITATIONS:START -->
## Scientific citations

These are the primary sources for the methods used in this workflow. Cite the relevant paper as well as JLinAlg when reporting results.

- [Geir Kjetil Sandve et al. (2013) — Ten Simple Rules for Reproducible Computational Research](../CITATIONS.md#sandve-reproducibility-2013)
- [W3C Provenance Working Group (2013) — PROV-O: The PROV Ontology](../CITATIONS.md#w3c-prov-o-2013)

[Search the complete scientific bibliography](https://robbyjo.github.io/JLinAlg/citations.html).
<!-- SCIENTIFIC-CITATIONS:END -->
