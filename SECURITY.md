# Security Policy

## Reporting a vulnerability

Please report vulnerabilities privately via GitHub's security advisories:
[**Report a vulnerability**](https://github.com/rpg32/caldyr/security/advisories/new)
(Security tab → "Report a vulnerability"). Please do not open a public issue
for anything you believe is exploitable.

You can expect an acknowledgement within a few days. Caldyr is maintained by a
small team, so please allow a reasonable window for a fix before public
disclosure — we'll keep you updated in the advisory thread.

## Supported versions

Caldyr is pre-1.0: only the latest release (and `main`) receives security
fixes.

## Scope and threat model

Caldyr is a **local-first** tool: the API server is designed to bind to
`127.0.0.1` and serve a single local user. Reports we especially care about:

- Code execution or file access triggered by opening a shared `.flow` file
  (these are untrusted input — the UI and engine must treat them as such).
- Anything reachable by a malicious web page while the local API is running
  (cross-origin/WebSocket attacks against `localhost`).
- Vulnerabilities in the docs site or CI/release pipeline.

Denial-of-service against your own local instance (e.g. asking your own
machine for an enormous simulation) is generally out of scope. Exposing the
API on a network interface is not a supported configuration; hardening for a
hosted multi-tenant deployment is tracked separately and is not covered by
this policy yet.
