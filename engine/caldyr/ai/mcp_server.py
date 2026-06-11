"""Expose the Caldyr tools as an MCP server (stdio).

This is how MCP clients — Codex CLI, Claude Desktop, etc. — use Caldyr natively:
they connect to this server and call the same typed tools the in-process agent
uses, against one shared :class:`AgentSession`. No LLM lives here; the client is
the model.

    # register with Codex CLI:
    codex mcp add caldyr -- python -m caldyr.ai.mcp_server
    # then ask Codex to "build an ammonia loop with Caldyr and cost it".

Run standalone for a smoke check: ``python -m caldyr.ai.mcp_server``.
"""
from __future__ import annotations

import json

import mcp.types as types
from mcp.server.lowlevel import Server
from mcp.server.stdio import stdio_server

from .session import AgentSession
from .tools import TOOLS, dispatch


def build_server() -> tuple[Server, AgentSession]:
    """Construct the MCP server and the session its tools mutate."""
    server: Server = Server("caldyr")
    session = AgentSession()

    @server.list_tools()
    async def list_tools() -> list[types.Tool]:
        return [types.Tool(name=t.name, description=t.description, inputSchema=t.input_schema)
                for t in TOOLS]

    @server.call_tool()
    async def call_tool(name: str, arguments: dict) -> list[types.TextContent]:
        out = dispatch(session, name, arguments or {})
        return [types.TextContent(type="text", text=json.dumps(out))]

    return server, session


async def _main() -> None:
    server, _ = build_server()
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())


def main() -> None:
    import anyio
    anyio.run(_main)


if __name__ == "__main__":
    main()
