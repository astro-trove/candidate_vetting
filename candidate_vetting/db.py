"""
Util function(s) for queries against the catalog database
"""
from contextlib import contextmanager
from django.db import connections, router, transaction

@contextmanager
def jit_disabled(model):
    """Disable PostgreSQL JIT for queries against `model`'s database.

    The alias is resolved through the project's DATABASE_ROUTERS
    (CatalogRouter), so this works in any project that uses this app.
    SET LOCAL reverts automatically when the transaction ends.
    """
    alias = router.db_for_read(model)
    conn = connections[alias]
    with transaction.atomic(using=alias):
        if conn.vendor == "postgresql":  # no-op on e.g. SQLite test DBs
            with conn.cursor() as cur:
                cur.execute("SET LOCAL jit = off")
        yield
