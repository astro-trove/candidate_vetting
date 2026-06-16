class CatalogRouter:
    def db_for_read(self, model, **hints):
        if model._meta.app_label == 'candidate_vetting':
            return 'catalogs'
        return None

    def db_for_write(self, model, **hints):
        if model._meta.app_label == 'candidate_vetting':
            return 'catalogs'
        return None

    def allow_migrate(self, db, app_label, model_name=None, **hints):
        if app_label == 'candidate_vetting':
            return db == 'catalogs'
        return db == 'default'
