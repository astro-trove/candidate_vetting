"""
Util function(s) for queries against the catalog database
"""
from contextlib import contextmanager
from django.db import connections, router

@contextmanager
def jit_disabled(model):
    """Disable PostgreSQL JIT for queries against `model`'s database.

    The alias is resolved through the project's DATABASE_ROUTERS
    (CatalogRouter), so this works in any project that uses this app.

    JIT is turned off with a session-level SET, issued once per database
    connection (Django replaces `conn.connection` on reconnect, so a new
    connection gets the SET again). This avoids the extra round trips of a
    per-query transaction + SET LOCAL, which dominate the cost of the simple
    q3c-indexed catalog queries on a remote database. JIT therefore stays off
    for every query on the catalog connection after the first use.
    """
    alias = router.db_for_read(model)
    conn = connections[alias]
    if conn.vendor == "postgresql":  # no-op on e.g. SQLite test DBs
        conn.ensure_connection()
        if getattr(conn, "_jit_off_for", None) is not conn.connection:
            with conn.cursor() as cur:
                cur.execute("SET jit = off")
            # inside a transaction the SET could still be rolled back, so only
            # remember it when it has taken effect immediately (autocommit)
            if not conn.in_atomic_block:
                conn._jit_off_for = conn.connection
    yield
